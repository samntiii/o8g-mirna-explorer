#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import datasets
from transformers import AutoTokenizer, AutoModelForSequenceClassification, Trainer, TrainingArguments
import multimolecule  # register rnabert
import torch

# 1. Load Data
dataset = datasets.load_dataset('csv', data_files={'train': 'train_0.9.csv', 'validation': 'valid_0.9.csv'})

# 2. Replace T with U
dataset = dataset.map(lambda e: {'sequence': [seq.replace('T', 'U') for seq in e['sequence']]}, batched=True)

# 3. Load tokenizer and models
model_name = 'multimolecule/rnabert'
tokenizer = AutoTokenizer.from_pretrained(model_name, bos_token=None, eos_token=None)
model = AutoModelForSequenceClassification.from_pretrained(model_name, num_labels=2)

# 4. Data tokenization
def preprocess(batch):
    encoding = tokenizer(batch['sequence'], padding='max_length', truncation=True, max_length=256)
    
    # 把标签转为 one-hot，例如 1 -> [0,1], 0 -> [1,0]
    labels = batch['label']
    one_hot_labels = [[1,0] if l == 0 else [0,1] for l in labels]
    
    encoding['labels'] = one_hot_labels
    return encoding

dataset = dataset.map(preprocess, batched=True)
dataset.set_format(type='torch', columns=['input_ids', 'attention_mask', 'labels'])

# 5. Train
training_args = TrainingArguments(
    output_dir='./rnabert-finetuned_0.8',
    num_train_epochs=5,
    per_device_train_batch_size=16,
    per_device_eval_batch_size=16,
    evaluation_strategy='epoch',
    save_strategy='epoch',
    logging_steps=50,
    learning_rate=2e-5,
    load_best_model_at_end=True,
    metric_for_best_model='accuracy',
)

# 6. Evaluation
import numpy as np
from sklearn.metrics import accuracy_score, f1_score, precision_score, recall_score, roc_auc_score

def compute_metrics(pred):
    # 模型输出（logits） → 预测类别索引
    preds = np.argmax(pred.predictions, axis=1)
    
    # One-hot → 真实标签索引
    labels = np.argmax(pred.label_ids, axis=1)
    
    # AUC 需要概率值（我们这里取 pred.probs[:, 1]）
    try:
        probs = pred.predictions[:, 1]  # 第二列是预测为1的概率
        auc = roc_auc_score(labels, probs)
    except:
        auc = -1  # 防止小batch出错

    return {
        'accuracy': accuracy_score(labels, preds),
        'precision': precision_score(labels, preds),
        'recall': recall_score(labels, preds),   # sensitivity
        'f1': f1_score(labels, preds),
        'auc': auc
    }

# 7. Trainer
trainer = Trainer(
    model=model,
    args=training_args,
    train_dataset=dataset['train'],
    eval_dataset=dataset['validation'],
    compute_metrics=compute_metrics,
)

# 8. Start training and save
trainer.train()
trainer.evaluate()

tokenizer.save_pretrained('./rnabert-finetuned_0.9')

trainer.model.save_pretrained('./rnabert-finetuned_0.9')

