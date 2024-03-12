import sys
sys.path.append('scripts/smiles-gpt/')

import smiles_gpt as gpt
import pandas as pd

# load pre-trained smiles-gpt models and tokenizer
from transformers import GPT2Config, GPT2LMHeadModel, PreTrainedTokenizerFastche 
checkpoint = "checkpoints/benchmark-5m"
config = GPT2Config.from_pretrained(checkpoint)
model = GPT2LMHeadModel.from_pretrained(checkpoint)
tokenizer = PreTrainedTokenizerFast.from_pretrained(checkpoint)

# add adapter
model.add_adapter("MAPT_adapter")
model.train_adapter("MAPT_adapter")
model.set_active_adapters("MAPT_adapter") # freeze all adapters except this one

# call the autoregressive causal language model
batch_size = 256

autoreg_model = gpt.GPT2LitModel(
    model, 
    batch_size=batch_size,
    learning_rate=1e-4,
    weight_decay=0.01,
    adam_eps=1e-8,
    adam_betas=(0.9, 0.999),
    scheduler_T_max=150_000)


# prepare training data
data_path = './data_preprocessing/actives_list.csv'
dataset = gpt.LMDataModule(data_path, tokenizer,
                           batch_size=batch_size,
                           num_workers=32)

# model training
from pytorch_lightning import Trainer

trainer = Trainer(max_epochs=50, gpus=4)
trainer.fit(autoreg_model, dataset)

autoreg_model.save_pretrained('/trained_models')