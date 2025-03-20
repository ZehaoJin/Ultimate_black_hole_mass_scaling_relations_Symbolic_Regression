import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import re
import scipy.optimize as opt
from tqdm import tqdm
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr
from sympy import lambdify
import torch
from torch import nn, optim
from torch.autograd.functional  import hessian 

from torch import log10
from torch import exp
from torch import log

from scipy import stats
from rpy2.robjects import r, FloatVector


r('library(goftest)')

number_matching_pattern = r"(?<![a-zA-Z0-9_.])[+-]?(\d+\.\d+|\.\d+|\d+\.|\d+)(?:[eE][-+]?\d+)?"

def simplify_equation(equation):
    raw_equation = equation
    # replace the power sign
    raw_equation = raw_equation.replace('^','**')
    # evaluate the exp()
    equation = parse_expr(raw_equation).evalf()
    # swap log10(constant) with log(constant)/log(10)
    equation = re.sub(r'log10\((.*?)\)', r'log(\1)/log(10.0)', str(equation))
    #print(equation)
    # evaluate the log()
    equation = str(parse_expr(equation).evalf())
    return equation

def str2equ(equation):
    return lambdify(list(dict.fromkeys(re.findall(r'\bx\d+',equation))),equation)


r('library(goftest)')

number_matching_pattern = r"(?<![a-zA-Z0-9_.])[+-]?(\d+\.\d+|\.\d+|\d+\.|\d+)(?:[eE][-+]?\d+)?"

def general_fit(equation,df,verbose=True, simplification=True, intrinsic_scatter_guess=None,
                initial_lr=0.05, max_epoch = 10000, 
                lr_decay='linear', linear_decay_start = None,
                plateau_decay_factor=0.5,plateau_patience=10, plateau_threshold=1e-6):
    df=df.copy()
    
    ### simplify the equation
    if simplification:
        equation = simplify_equation(equation)
        if verbose:
            print('Initial simplification:',equation)
    initial_equation = equation

    ### Get y and y_std
    # check if df contains y
    if 'y' not in df.columns:
        raise ValueError('DataFrame must contain columns named "y".')
    # if no y_std, create a column of zeros
    if 'y_std' not in df.columns:
        print('Warning: DataFrame does not contain a column named "y_std", will assume zero error.')
        df['y_std'] = 0
    # Y to tensor
    Y = torch.tensor(df['y'].values, dtype=torch.float64)
    Y_std = torch.tensor(df['y_std'].values, dtype=torch.float64)


    ### Get x index and convert to fitting format
    # Find all unique matches first to avoid duplicates and incorrect indexing
    constants = list(set(re.findall(number_matching_pattern, equation)))
    # Sort constants by their length in descending order to replace longer numbers first, preventing partial replacement issues
    constants.sort(key=len, reverse=True)

    def replace_with_variable(match):
        # Find the matched number in the constants list and get its index
        number = match.group(1)
        index = constants.index(number)

        # Check if the match includes a minus sign
        if match.group(0).startswith('-'):
            sign = '-'
        else:
            sign = ''
        return f'{sign}self.p{index}'
    # replace the number with p0, p1, p2, ...
    equation = re.sub(number_matching_pattern, replace_with_variable, equation)
    if verbose:
        print('Fitting format:',equation)

    # get x_indexs
    x_indexs = re.findall(r'x(\d+)', equation)
    x_indexs = list(dict.fromkeys(x_indexs))
    x_indexs = [int(x_index) for x_index in x_indexs]
    if verbose:
        print('x_indexs:',x_indexs)

    # check if contains x
    for x_index in x_indexs:
        if 'x'+str(x_index) not in df.columns:
            raise ValueError(f'DataFrame must contain columns named "x{x_index}".')
    # if no x_std, create a column of zeros
    for x_index in x_indexs:
        if 'x'+str(x_index)+'_std' not in df.columns:
            print(f'Warning: DataFrame does not contain a column named "x{x_index}_std", will assume zero error.')
            df['x'+str(x_index)+'_std'] = 0

    # construct tensors for x and x_std
    x=[]
    x_stds=[]
    for x_index in x_indexs:
        x.append(torch.tensor(df['x'+str(x_index)].values, dtype=torch.float64))
        x_stds.append(torch.tensor(df['x'+str(x_index)+'_std'].values, dtype=torch.float64))

    ### Get the data size, number of unique variables and number of constants
    data_size=len(df)
    unique_number_variables=len(x_indexs)
    number_constants=len(constants)
    # if intrinsic_scatter_guess is not given, set it to initial RMSE error of the prediction
    if intrinsic_scatter_guess == None:
        inital_equation_func = str2equ(initial_equation)
        intrinsic_scatter_guess = torch.sqrt(torch.mean((Y - inital_equation_func(*x))**2)).detach().numpy()
        
    if verbose:
        print('data_size:',data_size)
        print('unique_number_variables:',unique_number_variables)
        print('number_constants:',number_constants)
        print('intrinsic_scatter_guess:',intrinsic_scatter_guess)
    


    ### Define the model of the equation
    class Model(nn.Module):
        
        def __init__(self):
            # assign initial values to parameters self.p[i], self.x[i] and self.intrinsic_scatter
            super(Model, self).__init__()
            # get self.p0, self.p1, self.p2, ...
            for i in range(number_constants):
                setattr(self, f'p{i}', nn.Parameter(torch.tensor(float(constants[i]), dtype=torch.float64)))
            # get self.x0, self.x1, self.x2, ...
            for i,x_index in enumerate(x_indexs):
                # fittable x if x has error
                if x_stds[i].sum() != 0:
                    setattr(self, f'x{x_index}', nn.Parameter(x[i].clone()))
                # non-fittable x if x has no error
                else:
                    setattr(self, f'x{x_index}', x[i].clone())
            # self.intrinsic scatter
            self.intrinsic_scatter = nn.Parameter(torch.tensor(intrinsic_scatter_guess, dtype=torch.float64))
        
        def forward(self,i):
            # calculate the model prediction
            # replace x[i] with self.x[i]
            eq = equation
            for x_index in x_indexs:
                eq = re.sub(r'\bx'+str(x_index)+r'\b', f'self.x{x_index}[i]', eq)
            
            # evaluate the equation
            return eval(eq)
        
    model = Model()
        
    ### Define the loss function
    def loglikelihood():
        # term0
        term0 = torch.log(torch.tensor(2 * torch.pi)) * data_size * (unique_number_variables + 1)

        #term1
        term1 = (torch.log(Y_std**2 + model.intrinsic_scatter**2)).sum()
        for j in range(len(x_indexs)):
            if x_stds[j].sum() != 0:
                term1 += (torch.log(x_stds[j]**2)).sum()

        #term2
        term2 = ((Y - model(torch.arange(data_size)))**2 / (Y_std**2 + model.intrinsic_scatter**2)).sum()
        
        #term3
        term3 = 0
        for j in range(len(x_indexs)):
            if x_stds[j].sum() != 0:
                term3 += (((x[j] - getattr(model, f'x{x_indexs[j]}')) / x_stds[j])**2).sum()

        return term0 + term1 + term2 + term3

    
    ### Define the optimizer
    optimizer = optim.Adam(model.parameters(), lr=initial_lr)

    ### lr decay
    if lr_decay == 'linear':
        # default decay starting from max_epoch/2
        if linear_decay_start == None:
            linear_decay_start = max_epoch/2
        # linear decay function
        def lr_lambda(epoch):
            total_decay_epochs = max_epoch - linear_decay_start
            if epoch < linear_decay_start:
                return 1.0
            else:
                decay_ratio = (epoch - linear_decay_start) / total_decay_epochs 
                return max(0.0, 1.0 - decay_ratio)
    
        scheduler_linear = optim.lr_scheduler.LambdaLR(optimizer, lr_lambda=lr_lambda)


    elif lr_decay == 'plateau':
        scheduler_plateau = optim.lr_scheduler.ReduceLROnPlateau(optimizer,  factor=plateau_decay_factor, patience=plateau_patience,threshold=plateau_threshold)
    else:
        raise ValueError('lr_decay must be "linear" or "plateau".')
    
   
    ### Training loop
    for epoch in range(max_epoch):
        optimizer.zero_grad()
        loss = loglikelihood()
        if epoch == 0:
            initial_loss = loss.item()
        loss.backward()
        optimizer.step()
        if lr_decay == 'linear':
            scheduler_linear.step()
        else:
            scheduler_plateau.step(loss)
        
        if epoch % 100 == 0 and verbose:
            current_lr = optimizer.param_groups[0]['lr']
            print(f'Epoch {epoch}, LR: {current_lr:.5f}, Loss: {loss.item()}', 'p:', [getattr(model, f'p{i}').item() for i in range(len(constants))],'intrinsic scatter:' ,model.intrinsic_scatter.item())
    
    # check if the loss is converged
    if loss.item() > initial_loss or torch.isnan(loss):
        print('Warning! Loss did not converge. Please try again with different learning rate.')
        raise ValueError('Loss did not converge. Please try again with different learning rate.')

    ### Get the 2nd derivative of the loss function at last epoch
    def compute_hessian(model):
        # Collect initial parameters to get shapes and sizes
        params = list(model.parameters())
        shapes = [p.shape for p in params]
        numels = [p.numel() for p in params]

        # Define a function that takes a flat tensor and returns the loss
        def flat_loss(flat_params):
            # Split the flat tensor into individual parameters
            split_params = []
            offset = 0
            for shape, numel in zip(shapes, numels):
                param = flat_params[offset:offset + numel].view(shape)
                split_params.append(param)
                offset += numel

            # Create a dictionary of parameters for dynamic access
            param_dict = {name: param for name, param in zip(model.state_dict(), split_params)}
            # also include non-fittable x
            for j,x_index in enumerate(x_indexs):
                if x_stds[j].sum() == 0:
                    param_dict[f'x{x_index}'] = x[j]


            eq = equation
            #print(eq)
            for x_index in x_indexs:
                eq = re.sub(r'\bx'+str(x_index)+r'\b', f'param_dict["x{x_index}"]', eq)
            for i_index in range(len(constants)):
                eq = re.sub(r'self.p'+str(i_index)+r'\b', f'param_dict["p{i_index}"]', eq)
            #print(eq)
            pred=eval(eq)
            #print(pred)

            # Manually compute the model's output using the split parameters
            # Re-implement the loss calculation using param_dict
            term0 = torch.log(torch.tensor(2 * torch.pi)) * data_size * (unique_number_variables + 1)
            
            term1 = (torch.log(Y_std**2 + param_dict['intrinsic_scatter']**2)).sum()
            for j in range(len(x_indexs)):
                if x_stds[j].sum() != 0:
                    term1 += (torch.log(x_stds[j]**2)).sum()

            term2 = ((Y - pred)**2 / (Y_std**2 + param_dict['intrinsic_scatter']**2)).sum()

            term3 = 0
            for j in range(len(x_indexs)):
                if x_stds[j].sum() != 0:
                    term3 += (((x[j] - param_dict[f'x{x_indexs[j]}']) / x_stds[j])**2).sum()
            
            #print(term0 + term1 + term2 + term3)
            return term0 + term1 + term2 + term3


        # Compute Hessian with respect to the flattened initial parameters
        initial_flat = torch.cat([p.detach().flatten() for p in params])
        hessian_matrix = hessian(flat_loss, initial_flat)
        
        return hessian_matrix.numpy()

    hessian_matrix = compute_hessian(model)/(2*data_size)

    # Get cholesky decomposition
    try:
        L = np.linalg.cholesky(hessian_matrix)
        # get log determinant of the hessian with cholesky decomposition
        logdet = 2 * np.sum(np.log(L.diagonal()))

        # inverse of the hessian matrix gives the covariance matrix
        P=np.linalg.inv(L).transpose()
        inv_hessian_matrix = np.dot(P,P.transpose())
        sqrt_cov=np.sqrt(np.diag(inv_hessian_matrix))

    except np.linalg.LinAlgError:
        print('Warning! cholesky decomposition has non-positive diagonal elements')
        logdet = np.nan
        sqrt_cov = np.nan

    # get the string of final fit equation
    optimized_equation = equation
    for i in range(len(constants)):
        optimized_equation = re.sub(r'self.p'+str(i)+r'\b', f'{getattr(model, f"p{i}").item()}', optimized_equation)
    if verbose:
        print('optimized equation:',optimized_equation)

      
    # Get the optimized parameters
    optimized_params = {name: param.data for name, param in model.named_parameters()}
    if verbose:
        print('optimized parameters:',optimized_params)

    # residual
    residual_array=((Y - model(torch.arange(data_size))) / torch.sqrt(Y_std**2 + model.intrinsic_scatter**2)).detach().numpy()
    
    dim = 1 # dim of Y
    for j in range(len(x_indexs)):
        if x_stds[j].sum() != 0:
            residual_array=np.vstack((residual_array,((x[j] - getattr(model, f'x{x_indexs[j]}')) / x_stds[j]).detach().numpy()))
            dim += 1 # dim +1 for each x with error
    
    # test the normality with Anderson of the residuals from a N(0,1) distribution
    if dim >1:
        p_value=np.zeros(residual_array.shape[0]+1)
    else:
        p_value=np.zeros(1)
    # full flattened residual
    r_data = FloatVector(residual_array.flatten())
    result = r['ad.test'](r_data,null='pnorm',mean=0,sd=1)
    p_value[-1] = result[1][0]
    # residuals for each parameter
    if dim >1:
        for j in range(residual_array.shape[0]):
            # Convert Python data to R vector
            r_data = FloatVector(residual_array[j])

            # Perform Anderson test (R's Anderson-Darling test)
            result = r['ad.test'](r_data,null='pnorm',mean=0,sd=1)
            p_value[j] = result[1][0]
        


    # BIC
    N_para=np.max([number_constants,unique_number_variables])+1
    for j in range(len(x_indexs)):
        if x_stds[j].sum() != 0:
            N_para += data_size
    if verbose:
        print('N_para: ',N_para, 'dim: ',dim)
    # BIC = loss.item() + N_para * np.log(data_size) - N_para * np.log(2 * np.pi)
    BIC = loss.item() + N_para * np.log(data_size)
    BIC_i = loss.item() + N_para * np.log(data_size) + np.log(logdet) - N_para * np.log(2 * np.pi)

    result = {}
    result['fitted_equation'] = optimized_equation
    result['logdet'] = logdet
    result['sqrt_cov'] = sqrt_cov
    result['loss'] = loss.item()
    result['optimized_params'] = optimized_params
    result['N_para'] = N_para
    result['residual_array'] = residual_array
    result['ad_p_value_array'] = p_value
    result['BIC'] = BIC
    result['BIC_i'] = BIC_i

    return result
    