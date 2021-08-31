import os, subprocess

for root, dirs, files in os.walk(".", topdown=False):
   for name in files:
      name = os.path.join(root, name)
      if name.endswith('.ipynb'):
          print(name, 'ok')
          cmd = f"jupyter nbconvert --clear-output --inplace {name}"
          subprocess.call(cmd, shell=True)
