import numpy as np
from sklearn.linear_model import LinearRegression
def estimate_parameters(X, y):
    """
    Estimate the parameters of a linear regression model using least squares.

    Parameters:
    X (numpy array): The input features, shape (n_samples, n_features).
    y (numpy array): The target values, shape (n_samples,).

    Returns:
    numpy array: The estimated parameters (coefficients), shape (n_features,).
    """
    # Create a LinearRegression model
    model = LinearRegression()
    
    # Fit the model to the data
    model.fit(X, y)
    
    # Return the estimated coefficients
    return model.coef_
