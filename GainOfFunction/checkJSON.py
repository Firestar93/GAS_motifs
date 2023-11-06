import json

def print_json_structure(d, indent=0):
    for key, value in d.items():
        print('  ' * indent + str(key))
        if isinstance(value, dict):
            print_json_structure(value, indent+1)
        elif isinstance(value, list) and value and isinstance(value[0], dict):
            print('  ' * (indent+1) + "List of dictionaries with keys:")
            print_json_structure(value[0], indent+2)

def analyze_json_file(file_path):
    with open(file_path, 'r') as file:
        data = json.load(file)
        print_json_structure(data)

# Replace 'path_to_your_json_file.json' with the path to your JSON file
analyze_json_file('C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs\\GainOfFunction\\stuff\\response_1699101962308.json')
