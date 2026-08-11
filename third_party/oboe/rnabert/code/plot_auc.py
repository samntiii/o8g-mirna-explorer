#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#!/usr/bin/env python
# coding: utf-8

import os
import datasets
import torch
from transformers import AutoTokenizer, AutoModelForSequenceClassification
from sklearn.metrics import roc_curve, auc
import matplotlib.pyplot as plt
import numpy as np

import multimolecule  

# 1. Load dataset
dataset = datasets.load_dataset('csv', data_files={'validation': 'valid_0.9.csv'})
dataset = dataset.map(lambda e: {'sequence': [seq.replace('T', 'U') for seq in e['sequence']]}, batched=True)

# 2. Load tokenizer and models
model_name_or_path = './rnabert-finetuned_0.9'  # 使用finetuned模型
tokenizer = AutoTokenizer.from_pretrained(model_name_or_path, bos_token=None, eos_token=None)
model = AutoModelForSequenceClassification.from_pretrained(model_name_or_path)
model.eval()  # 设置为eval模式

# 3. Tokenize
def preprocess(batch):
    encoding = tokenizer(batch['sequence'], padding='max_length', truncation=True, max_length=256)
    encoding['label'] = batch['label']
    return encoding

dataset = dataset.map(preprocess, batched=True)
dataset.set_format(type='torch', columns=['input_ids', 'attention_mask', 'label'])

# 4. Use model to predict
all_labels = []
all_probs = []

with torch.no_grad():
    for batch in torch.utils.data.DataLoader(dataset['validation'], batch_size=16):
        inputs = {k: v for k, v in batch.items() if k in ['input_ids', 'attention_mask']}
        outputs = model(**inputs)
        probs = torch.softmax(outputs.logits, dim=1)[:, 1]  # 取出预测为1的概率
        all_probs.extend(probs.cpu().numpy())
        all_labels.extend(batch['label'].cpu().numpy())

# 5. Plot ROC figures
fpr, tpr, thresholds = roc_curve(all_labels, all_probs)
roc_auc = auc(fpr, tpr)

# 5.1 save ROC files
model_tag = "RNABERT_0.9"
save_dir = "./roc_data"
os.makedirs(save_dir, exist_ok=True)  # 自动创建文件夹

np.save(f"{save_dir}/fpr_{model_tag}.npy", fpr)
np.save(f"{save_dir}/tpr_{model_tag}.npy", tpr)

with open(f"{save_dir}/auc_{model_tag}.txt", "w") as f:
    f.write(str(roc_auc))

print("保存成功！文件位置：", save_dir)

plt.figure(figsize=(8, 6))
plt.plot(fpr, tpr, color='darkorange', lw=2, label=f'ROC curve (AUC = {roc_auc:.4f})')
plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
plt.xlim([-0.05, 1.05])
plt.ylim([-0.05, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver Operating Characteristic')
plt.legend(loc="lower right")
plt.grid(True)
plt.savefig('roc_curve.png', dpi=300)
plt.show()

