import json

# File paths
input_file = '/shared/results/michalis_exosomes/smallRNA/target_prediction/iri_random_vs_hsa_UTR_112/iri_random_target_conservation.json'
output_file = '/shared/results/michalis_exosomes/smallRNA/target_prediction/iri_random_vs_hsa_UTR_112/iri_random_target_conservation_depurated.json'

# List of identifiers to remove
identifiers_to_remove = [
    "iri-random-3",
    "iri-random-7",
    "iri-random-100",
    "iri-random-2",
    "iri-random-12",
    "iri-random-15",
    "iri-random-16",
    "iri-random-22",
    "iri-random-24",
    "iri-random-25",
    "iri-random-31",
    "iri-random-41",
    "iri-random-46",
    "iri-random-51",
    "iri-random-56",
    "iri-random-65",
    "iri-random-67",
    "iri-random-80",
    "iri-random-91",
    "iri-random-95"
]

# Load the JSON file
with open(input_file, 'r') as file:
    data = json.load(file)

# Remove specified identifiers
for identifier in identifiers_to_remove:
    if identifier in data:
        del data[identifier]

# Write the updated data back to a new JSON file
with open(output_file, 'w') as file:
    json.dump(data, file, indent=4)  # Using indent for pretty printing

print("Depuration complete. Updated JSON saved to:", output_file)

