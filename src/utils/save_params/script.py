import re
import yaml
import base64

## VIASH START
par = {
    "id": "sample_one",
    "params_yaml": "cGFyYW1zX3lhbWw6IHt9Cg==",
    "output": "output.yaml"
}
## VIASH END

class Dumper(yaml.Dumper):
    def increase_indent(self, flow=False, indentless=False):
        return super(Dumper, self).increase_indent(flow, False)

def decode_params_yaml(encoded_yaml):
    # Step 1: Decode from Base64
    yaml_bytes = base64.b64decode(encoded_yaml)
    
    # Step 2: Convert bytes to string
    yaml_string = yaml_bytes.decode('utf-8')
    
    # Step 3: Extract pattern for Java path objects
    # Find pattern: !!sun.nio.fs.UnixPath /path/to/file
    pattern = r'!!sun\.nio\.fs\.UnixPath\s+([^\n]+)'
    
    # Replace with the actual path string (captured group)
    yaml_string = re.sub(pattern, r'\1', yaml_string)
    
    # Handle any remaining empty UnixPath objects
    yaml_string = yaml_string.replace('!!sun.nio.fs.UnixPath {}', '""')
    
    # Step 4: Parse YAML
    yaml_data = yaml.safe_load(yaml_string)
    
    return yaml_data

def replace_id(value, sample_id):
    if isinstance(value, str):
        return value.replace('$id', sample_id)
    elif isinstance(value, list):
        return [replace_id(item, sample_id) for item in value]
    return value

print(par['params_yaml'])

params = decode_params_yaml(par['params_yaml'])
for key, value in params.items():
    params[key] = replace_id(value, par["id"])

with open(par["output"], 'w') as f:
    yaml.dump(params, f, default_flow_style=False, Dumper=Dumper)

