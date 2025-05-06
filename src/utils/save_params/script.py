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

print(par['params_yaml'])

params = decode_params_yaml(par['params_yaml'])

with open(par["output"], 'w') as f:
    yaml.dump(params, f, default_flow_style=False, Dumper=Dumper)

