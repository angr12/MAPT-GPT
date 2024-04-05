# run hyperparam search
import optuna
import pytorch_lightning as pl

import smiles_gpt as gpt
import pandas as pd

from transformers import GPT2Config, GPT2LMHeadModel, PreTrainedTokenizerFast

from pytorch_lightning import Trainer
from pytorch_lightning.callbacks import ModelCheckpoint, EarlyStopping

def objective(trial):
    # hyperparams search space
    batch_size = trial.suggest_int("batch_size", 16, 256)
    learning_rate = trial.suggest_float("learning_rate", 1e-5, 1e-2)
    weight_decay = trial.suggest_float("weight_decay", 0.0, 0.1)
    
    # load pre-trained smiles-gpt and tokenizer
    checkpoint = "checkpoints/benchmark-5m"
    config = GPT2Config.from_pretrained(checkpoint)
    model = GPT2LMHeadModel.from_pretrained(checkpoint)
    tokenizer = PreTrainedTokenizerFast.from_pretrained(checkpoint)

    # Set pad token
    tokenizer.pad_token = "<pad>"
    
    # load training data
    data_path = 'data_preprocessing/actives_list.csv'
    dataset = gpt.LMDataModule(data_path, tokenizer,
                               batch_size=batch_size,
                               num_workers=32)

    # add adapter
    model.add_adapter("MAPT_hyperparam_search")
    model.train_adapter("MAPT_hyperparam_search")
    model.set_active_adapters("MAPT_hyperparam_search") # freeze all adapters except this one

    transformer = gpt.GPT2LitModel(
        transformer=model, 
        batch_size=batch_size,
        learning_rate=learning_rate,
        final_learning_rate=1e-5,
        weight_decay=weight_decay,
        adam_eps=1e-8,
        adam_betas=(0.9, 0.999),
        scheduler_T_max=150_000
    )
    
    trainer = Trainer(
        strategy="ddp",
        max_epochs=35,
        min_epochs=15,
        val_check_interval=0.4,
        limit_train_batches=0.5,
        log_every_n_steps=10,
        gpus=-1
    )
    
    trainer.fit(transformer, dataset)
    
    ppl_epoch = trainer.callback_metrics['ppl_epoch']
    print(f'ppl_epoch: {ppl_epoch}') # debug line
    return ppl_epoch

if __name__ == "__main__":
    study = optuna.create_study(direction="minimize") # minmise perplexity
    study.optimize(objective, n_trials=100)
    
    print("Number of finished trials: {}".format(len(study.trials)))
    
    print("Best trial:")
    trial = study.best_trial
    
    print("  Value: {}".format(trial.value))
    
    print(' Params: ')
    for key, value in trial.params.items():
        print("    {}: {}".format(key, value))