#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#!/usr/bin/env python
# -*- coding: utf-8 -*-

import torch
import numpy as np
from transformers import AutoTokenizer, AutoModelForSequenceClassification
from Bio import SeqIO
import argparse

# -----------------------
# 模型路径：使用你finetune好的checkpoint
# -----------------------
MODEL_PATH = "./rnabert-finetuned_0.9/checkpoint-990"

# -----------------------
# 加载Tokenizer和Model
# -----------------------
print("Loading model and tokenizer...")
tokenizer = AutoTokenizer.from_pretrained(MODEL_PATH)
model = AutoModelForSequenceClassification.from_pretrained(MODEL_PATH)
model.eval()

# -----------------------
# 预测函数
# -----------------------
def predict_o8g(sequence):
    sequence = sequence.replace("T", "U") 
    inputs = tokenizer(sequence, return_tensors="pt", padding="max_length", truncation=True, max_length=256)

    with torch.no_grad():
        outputs = model(**inputs)
        logits = outputs.logits
        probs = torch.softmax(logits, dim=1).numpy()[0]

    pred_class = int(np.argmax(probs))  # 0 or 1
    pred_prob = float(probs[1])         # o8G 的概率
    return pred_class, pred_prob

# -----------------------
# 从FASTA中读取多个序列
# -----------------------
def predict_from_fasta(fasta_path, output_path=None):
    results = []

    for record in SeqIO.parse(fasta_path, "fasta"):
        seq = str(record.seq)
        label, prob = predict_o8g(seq)
        result = {
            "id": record.id,
            "prediction": "YES" if label == 1 else "NO",
            "probability": prob
        }
        results.append(result)

    # 打印结果
    for r in results:
        print(f">{r['id']}\nPrediction: {r['prediction']} | Probability: {r['probability']:.4f}")

    # 可选：保存为TSV
    if output_path:
        with open(output_path, 'w') as f:
            f.write("id\tprediction\tprobability\n")
            for r in results:
                f.write(f"{r['id']}\t{r['prediction']}\t{r['probability']:.4f}\n")
        print(f"\nResults saved to {output_path}")

# -----------------------
# CLI 接口
# -----------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Predict o8G modification from fasta sequences.")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file path.")
    parser.add_argument("-o", "--output", help="Output TSV file path (optional).")
    args = parser.parse_args()

    predict_from_fasta(args.input, args.output)

