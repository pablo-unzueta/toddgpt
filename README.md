# Todd-GPT

A LLM-agent for Excited State Chemistry Simulations

## Setup

This project uses Docker for easy setup and deployment. Follow these steps to get started:

### Prerequisites

1. Install Docker: https://docs.docker.com/get-docker/

### Using the pre-built Docker image

1. Pull the Docker image:
   ```
   docker pull pablounzueta/toddgpt:latest
   ```

2. Run the Docker container:
   ```
   docker run -it --rm \
     -e OPENAI_API_KEY=$OPENAI_API_KEY \
     -v $(pwd):/app/workdir \
     pablounzueta/toddgpt:latest
   ```

### Building the Docker image locally

If you want to build the image yourself:

1. Clone this repository:
   ```
   git clone https://github.com/pablo-unzueta/toddgpt.git
   cd toddgpt
   ```

2. Build the Docker image:
   ```
   docker build -t toddgpt:latest .
   ```

3. Run the Docker container:
   ```
   docker run -it --rm \
     -e OPENAI_API_KEY=$OPENAI_API_KEY \
     -v $(pwd):/app/workdir \
     toddgpt:latest
   ```

## Usage

1. When you run the Docker container, you'll be prompted to enter a question.
2. If your question involves a specific geometry file (e.g., an .xyz file), make sure it's in your current directory before running the Docker container.
3. The agent will process your question and provide a response based on the available information and files.

Example question:
```
What are the distances between atoms in a water molecule in h2o.xyz?
```

## Development

If you want to modify the code or contribute to the project, you can clone the repository and work on it locally. The project uses Poetry for dependency management, so you'll need to install Poetry and set up the virtual environment before making changes.

1. Clone the repository:
   ```
   git clone https://github.com/pablo-unzueta/toddgpt.git
   cd toddgpt
   ```

2. Install Poetry:
   ```
   pip install poetry
   ```

3. Set up the virtual environment:
   ```
   poetry shell
   ```

4. Install dependencies:
   ```
   poetry install
   ```

5. Make your changes to the code.

6. Build the Docker image:
   ```
   docker build -t toddgpt:latest .
   ```

7. Run the Docker container:
   ```
   docker run -it --rm \
     -e OPENAI_API_KEY=$OPENAI_API_KEY \
     -v $(pwd):/app/workdir \
     toddgpt:latest
   ```
