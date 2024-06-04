import smiles_gpt as gpt
import pandas as pd

# load pre-trained smiles-gpt models and tokenizer
from transformers import GPT2Config, GPT2LMHeadModel, PreTrainedTokenizerFast
checkpoint = "checkpoints/benchmark-5m"
config = GPT2Config.from_pretrained(checkpoint)
model = GPT2LMHeadModel.from_pretrained(checkpoint)
tokenizer = PreTrainedTokenizerFast.from_pretrained(checkpoint)

# Set pad token
tokenizer.pad_token = "<pad>"

# add adapter
model.add_adapter("MAPT_adapter")
model.train_adapter("MAPT_adapter")
model.set_active_adapters("MAPT_adapter") # freeze all adapters except this one

# call the autoregressive causal language model
batch_size = 128
 
autoreg_model = gpt.GPT2LitModel(
    transformer=model, 
    batch_size=batch_size,
    learning_rate=0.005,
    final_learning_rate=1e-5,
    weight_decay=4.6e-5,
    adam_eps=1e-8,
    adam_betas=(0.9, 0.999),
    scheduler_T_max=150_000)


# prepare training data
data_path = 'data_preprocessing/actives_list.csv'
dataset = gpt.LMDataModule(data_path, tokenizer,
                           batch_size=batch_size,
                           num_workers=32)

# model training
from pytorch_lightning import Trainer
from pytorch_lightning.callbacks import ModelCheckpoint, EarlyStopping

checkpoint_cb = ModelCheckpoint(dirpath='checkpoints/trained_models')

early_stopping = EarlyStopping(
    monitor='ppl_epoch',
    min_delta=0,
    patience = 3,
    verbose=True,
    mode='min'
)

trainer = Trainer(
    strategy="ddp",
    callbacks=[checkpoint_cb, early_stopping],
    max_epochs=35,
    min_epochs=15,
    val_check_interval=0.4,
    limit_train_batches=0.5,
    log_every_n_steps=10,
    gpus=-1
)
trainer.fit(autoreg_model, dataset)

autoreg_model.transformer.save_pretrained('checkpoints/trained_models/model/')