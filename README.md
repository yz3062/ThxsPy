# ThxsPy
This software calculates thorium excess, a constant flux proxy for ocean sediments, from ICP mass spectrometer results of uranium and thorium isotopes. It has a sister project, [PaxsPy](https://github.com/yz3062/PaxsPy), that calculates protactinium excess.

## Getting started

### Prerequisites

If you're new to python, I recommend you simply download [Anaconda](https://www.anaconda.com/download/), which includes all the prerequisites. Alternatively, download the followings individually:

Python 2

[scipy](https://www.scipy.org/)

[numpy](http://www.numpy.org/)

[pandas](https://pandas.pydata.org/)

### Installing

Click "Clone or download" on top right -> "Download ZIP". Unzip all files.

## Running the tests

Open Command Prompt in Windows or Terminal in MacOS, change to the directory where ThxsPy exsits, and type
```
python ThxsPy.py
```
You'll be prompted to confirm the version of Th spikes you're using

![alt text](/README_screenshots/spike_prompt_ThxsPy.JPG)

Click "Yes", and you'll then be prompted to select data files as well as a sample info file

![alt text](/README_screenshots/data_selection_ThxsPy.JPG)

Where you should double click "data" folder and select all the files in that folder

![alt text](/README_screenshots/data_select_all_ThxsPy.JPG)

Notice that alongside the data files, there's also a "sample_info.xlsx" file that looks like this

![alt text](/README_screenshots/sample_info_screenshot.JPG)

And Voila! Calculation is done and you're asked to save the output file

![alt text](/README_screenshots/save_output_ThxyPy.JPG)


## License

[BSD](https://opensource.org/licenses/BSD-2-Clause)
