import torch

def get_gpu_count():
    return torch.cuda.device_count()

print(f"Number of available GPUs: {get_gpu_count()}")