import smiles_gpt as gpt
from transformers import pipeline, GPT2Config, GPT2LMHeadModel, PreTrainedTokenizerFast
import pandas as pd

# load trained model
model = GPT2LMHeadModel.from_pretrained("checkpoints/trained_models/model")
tokenizer = PreTrainedTokenizerFast.from_pretrained('checkpoints/benchmark-5m')
model.set_active_adapters("MAPT_adapter")

print(f'Loaded adapter: {model.active_adapters}') # debug line

# Set pad token
tokenizer.pad_token = "<pad>"

# create molecule generation pipeline 
generator = pipeline(
    'text-generation', 
    model=model, 
    tokenizer=tokenizer, 
    do_sample=True,
    top_k=50,
    top_p=0.96,
    num_return_sequences=200
    )
output = generator('', max_length=50)

# for testing/debug
# for i in range(len(output)):
#     print(f'Generated molecule: {output[i]}')

# save generated molecules to csv
df = pd.DataFrame(output)
df.to_csv('generated_molecules/top_k.csv', index=False)