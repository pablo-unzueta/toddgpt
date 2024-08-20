import os
from toddgpt.agent import Agent

def main():
    print("Starting ToddGPT...")
    
    # Initialize the Agent with OpenAI API key
    api_key = os.environ.get("OPENAI_API_KEY")
    if not api_key:
        print("OPENAI_API_KEY not found in environment variables.")
        raise ValueError("Please set the OPENAI_API_KEY environment variable.")
    
    print("Initializing Agent...")
    agent = Agent("openai", api_key)
    executor = agent.get_executor()

    # Get the question from the user
    conversation = input("Please enter your question: ")

    # Check if the question mentions a file
    file_path = None
    words = conversation.split()
    for i, word in enumerate(words):
        if word.endswith('.xyz'):
            file_path = word
            # Convert the file path to the Docker container path
            docker_file_path = os.path.join('/app/workdir', os.path.basename(file_path))
            words[i] = docker_file_path
            break

    # Reconstruct the conversation with the updated file path
    conversation = ' '.join(words)

    # If a file was mentioned, check if it exists
    if file_path:
        if not os.path.exists(docker_file_path):
            print(f"Error: The file {file_path} does not exist in the current directory.")
            return

    print(f"Processing conversation: {conversation}")

    # Run the agent
    print("Running agent...")
    result = executor.invoke({"conversation": conversation})

    # Print the result
    print("Agent's response:")
    print(result['output'])

if __name__ == "__main__":
    main()