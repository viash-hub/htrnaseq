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
    yaml_bytes = base64.b64decode(encoded_yaml)
    yaml_string = yaml_bytes.decode('utf-8')
    yaml_data = yaml.safe_load(yaml_string)
    
    return yaml_data

# def replace_id(value, sample_id):
#     if isinstance(value, str):
#         return value.replace('$id', sample_id)
#     elif isinstance(value, list):
#         return [replace_id(item, sample_id) for item in value]
#     return value

print(par['params_yaml'])

params = decode_params_yaml(par['params_yaml'])
# for key, value in params.items():
#     params[key] = replace_id(value, par["id"])

with open(par["output"], 'w') as f:
    yaml.dump(params, f, default_flow_style=False, Dumper=Dumper)

