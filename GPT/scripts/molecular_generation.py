import smiles_gpt as gpt
from transformers import pipeline, GPT2Config, GPT2LMHeadModel, PreTrainedTokenizerFast

# import torch
# torch_device = "cuda" if torch.cuda.is_available() else "cpu"

# load trained model
model = GPT2LMHeadModel.from_pretrained("checkpoints/trained_models")
tokenizer = PreTrainedTokenizerFast.from_pretrained('checkpoints/benchmark-5m')
model.set_active_adapters("MAPT_adapter")

print(f'Loaded adapter: {model.active_adapters}') # debug line

# Set pad token
tokenizer.pad_token = "<pad>"

# create molecule generation pipeline 
generator = pipeline('text-generation', model=model, tokenizer=tokenizer)
output = generator('C', max_length=50)

print(output)