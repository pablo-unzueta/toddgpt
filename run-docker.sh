docker run -it --rm \
  -e OPENAI_API_KEY=$OPENAI_API_KEY \
  -e CHEMCLOUD_USER=$CHEMCLOUD_USER \
  -v $(pwd):/app/workdir \
  toddgpt:latest \
  python main.py