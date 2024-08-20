FROM python:3.9.16-slim

# Set working directory
WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    gcc \
    && rm -rf /var/lib/apt/lists/*

# Install poetry
RUN pip install poetry

# Copy poetry files
COPY pyproject.toml poetry.lock ./

# Install dependencies
RUN poetry config virtualenvs.create false \
  && poetry install --no-interaction --no-ansi

# Copy your application code
COPY . .

# Set Python path
ENV PYTHONPATH=/app/src

# Command to run your application
CMD ["python", "-m", "toddgpt"]