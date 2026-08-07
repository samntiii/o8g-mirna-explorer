# o8G-miRNA Retargeting Explorer — container image
FROM python:3.13-slim

WORKDIR /app

# system deps kept minimal; wheels cover numpy/scipy/pyarrow
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# app code + precomputed data
COPY . .

# unpack gene-sets if only the tarball shipped
RUN if [ ! -d genesets ] && [ -f genesets.tar.gz ]; then tar -xzf genesets.tar.gz; fi

# managed platforms inject $PORT; default to 8501 locally
ENV PORT=8501
EXPOSE 8501

# 0.0.0.0 so the container is reachable; headless off the browser auto-open
CMD streamlit run app.py \
    --server.port=${PORT} \
    --server.address=0.0.0.0 \
    --server.headless=true \
    --server.enableCORS=false \
    --server.enableXsrfProtection=false
