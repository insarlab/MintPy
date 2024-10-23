# This is a new file

import torch

def matrix_multiply(matrix_a, matrix_b):
    """
    Multiplies two matrices using PyTorch.

    Args:
    - matrix_a (list of list of floats): The first matrix.
    - matrix_b (list of list of floats): The second matrix.

    Returns:
    - result (torch.Tensor): The result of the matrix multiplication.
    """
    tensor_a = torch.tensor(matrix_a)
    tensor_b = torch.tensor(matrix_b)
    result = torch.matmul(tensor_a, tensor_b)
    return result