#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import torch
from torch.optim import AdamW
from torch.nn import functional as F
from transformers import AutoModelForSequenceClassification, get_linear_schedule_with_warmup
from preprocess import load_dataloader

def train_model(dataloader, model_name='zhihan1996/DNA_bert_6', num_labels=2, num_epochs=3, lr=2e-5):
    print("🚀 开始加载模型并进行训练...")
    
    device = torch.device("mps" if torch.backends.mps.is_available() else "cuda" if torch.cuda.is_available() else "cpu")
    model = AutoModelForSequenceClassification.from_pretrained(model_name, num_labels=num_labels)
    model.to(device)

    optimizer = AdamW(model.parameters(), lr=lr)
    total_steps = len(dataloader) * num_epochs
    scheduler = get_linear_schedule_with_warmup(optimizer, num_warmup_steps=0, num_training_steps=total_steps)

    for epoch in range(num_epochs):
        model.train()
        running_loss = 0.0
        for batch in dataloader:
            input_ids = batch['input_ids'].to(device)
            attention_mask = batch['attention_mask'].to(device)
            labels = batch['labels'].to(device)

            optimizer.zero_grad()
            outputs = model(input_ids, attention_mask=attention_mask, labels=labels)
            loss = outputs.loss
            loss.backward()
            optimizer.step()
            scheduler.step()

            running_loss += loss.item()

        avg_loss = running_loss / len(dataloader)
        print(f"📘 Epoch {epoch + 1}/{num_epochs} - Loss: {avg_loss:.4f}")

    save_path = './finetuned_dnabert_6mer_1'
    model.save_pretrained(save_path)
    print(f"✅ 模型保存成功：{save_path}")

if __name__ == "__main__":
    dataloader = load_dataloader()
    train_model(dataloader)

