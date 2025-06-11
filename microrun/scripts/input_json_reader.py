#! /usr/bin/env python3

import json

def read_json_file(file_path):
    with open(file_path, 'r') as file:
        data = json.load(file)
    for key, value in data.items():
        # Ensure values with spaces are enclosed in quotes
        print(f'{key}="{value}"')

if __name__ == "__main__":
    file_path = "input.json"
    read_json_file(file_path)