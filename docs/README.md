# Introduciton
SCAWP is a web app for easy-to-use single cell analysis. It aims to provide quick preprocessing, 
beatiful visualizations, and rich interactions for researcher without much experience in 
programming. 

Currently, SCAWP is under early development. This documentation will target at the potential 
developers of this project.   

SCAWP is a browser app server that uses Django as backend and Browser frontend. It can be 
deployed locally or remotely as Jupyter and accessed from the frontend. The app should be able to 
distribute as a Python package or a software. 
 
# Workflow for (ideal) the final product   
| Client | server |
| --- | --- |
| upload dataset in various forms to the server | read using the provided reader calls and save them as `h5ad` |
| building the preprocessing steps | parse the restAPI, resolve all the errors, and execute each step as one python instruction. Then save the results in the corresponding directory | 
| view, save, and export the results | serve all the results, including processed dataset, result report, plots | 
| save the preprocessing steps and apply on other datasets | save the processing calls and make adjustment for other datasets |
| view and interact with the data, including subseting and further computations | serve the plots and respond to any necessary computations requested |
  

# Tech Stacks
__server__ `Django`   
__processing/computation__ `scanpy` (core), other packages are also used but should be optional.  
__visualization__ `plotly`  
__frontend__ JS, HTML, CSS, the CSS template can be found on https://dashboardpack.com/theme-details/architectui-html-dashboard-free/.  
There are additional JS/CSS libraries but currently we use them directly from CDN for future changes.  


# Project Structure
The application is divided into several Django apps that corresponds to the resources and API calls for
different parts of the frontend functionality. 
 - `dataset/` uploading (saving), viewing, computing, exporting on the datasets
 - `process/` parse the restAPI calls, validate request,  and execute Python code. 
 Currently, we use __unsafe__ ways to achieve the maximum freedom. This app also functions 
 for all the management on the process history.
 - `plot/` handle the resources and computations for plot studio (real-time interactive plots) and 
 plots (static, saved plots).
 - `settings/` the user settings for supported function calls and plots.  
 
 # `iplot`
 `iplot` is a standalone package for interactive plotting. It works as an additional library over 
 other package's plotting modules. Typically, it will in the form of `iplot.[package_name].[function_call]`
 and for each module, it should have similar interfaces as the original function call. 
 
 # RestAPI for preprocessing
 When the user construct the preprocessing as building blocks. The frontend is essentially construct 
 an array of objects and send back to the server. Therefore, the supported computations should be 
 specified as follows so that the server can correctly parse the object and execute wanted code. 
 
 All function calls will be provided as a JSON object as follows, write a parser or manually add them. If you add them manually, 
 you can click the top right corner inside the portal and in the installed method section there is the wizard for adding methods. 
 
  `name:string, package:string` The exact function name and package name of the function.

 `description:string` a help string that explain what this function does

 `params:Array[Object]` necessary parameters for the function call excludes the first parameter, params will be an array of parameters. Each parameter will be a json object in the format specified

 `prerequisites: Array[Object] (optional)` if this method need to run after some other methods. Each prerequisite will be a json object in the form of name and package.


##### Specifications for parameter
 `name: string` parameter's name

 `type: [ text | number | option | bool ]` The type of the input,
  - `text` if the function parameter is a `str` or `List(str)`
  - `number` if the function parameter is `int`, `float` or other numerical type
  - `option` if the function parameter is one of several values
  - `bool` if the function parameter is a bool

`annotation: string (optional)` explaination for this parameter

`isList: bool (optional, type=text)`  if `List(str)` is expected for parameter input. Then the string will be splited by `, `

`options: Array[string] (required if type=option)` the options.

`default: string (optinoal)` the default value of the parameter.

`required: bool (optional)` if the parameter is required


# Modules to Implement
- more functions within plot studio (ideally something similar to BioTuring https://bioturing.com/).
- Smart suggestions for processing and better prerequisite check. Many visualizations or interactions 
requires further computation and adding new data attributes onto annData. Want to integrate the
preprocessing module with the plot studio so that it can be fired within plot studio and save the 
proper results and record.   
- `procedure`, a.k.a. a wrapper for a block of function calls, so that the user can create "functions" 
from the frontend and reuse them. 
- efficiency, the current response speed is not ideal enough, find ways to optimize, considering 
WebSocket. 

# Optional Suggestions for a Larger Scale
Because of the poor scalibility of dynamic languages, if we want to scale up this project to 
adapt more possibilities. I'm thinking to divide the app into frontend resources, restAPI backend,
and computation modules.  
By lowering the coupling, this project can then be developed in a larger scale, and we can achieve 
higher frontend stability by introducing React/TypeScript + NodeJS for the frontend resources, and 
lower delay on the computations by replacing Django with fastAPI (around 10x response speed of Django).   