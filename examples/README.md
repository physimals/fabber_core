# Example model template

This example provides a complete template module for a model of your own. To adapt it to your own model you will
need to do, as a minimum, the following:

 - Replace `sine` in all files with the name of your model (including `CMakeLists.txt`)

 - Rename files to replace `sine` with your model name
 
 - Edit `<model_name>_models.cc` to return the name of your model and its constructor function
 
 - Edit `fwdmodel_<model>.cc` and `fwdmodel_<model>.h` to define your model!

