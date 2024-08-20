docker run -it --rm \
  -e OPENAI_API_KEY=$OPENAI_API_KEY \
  -v $(pwd):/app/workdir \
  toddgpt:latest \
  python main.py