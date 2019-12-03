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

Included in this program is a test dataset to get you familiar with how it works. Open Command Prompt in Windows or Terminal in MacOS, change to the directory where ThxsPy exsits, and type
```
python ThxsPy.py
```
You'll be prompted to confirm the version of Th spikes you're using

![alt text](/README_screenshots/spike_prompt_ThxsPy.png)

Hit "Enter", and you'll be asked whether you want to inspect the data in figures

![alt text](/README_screenshots/inspect_figures_prompt.png)

Hit "Enter". Figures of all the ICPMS counts will be saved in the same folder as the input files. You can check the figures to see if there's anything abnormal, e.g. a spike in counts or trailing in counts midway. Notice that the all isotopes are plotted on the same y-axis, meaning you'll mostly see the variations in major isotopes like 238U and 232Th. You'll then select data files as well as a sample info file. Notice that the file selector window sometimes doesn't pop up and is open in the background.

![alt text](/README_screenshots/data_selection_ThxsPy.JPG)

Where you should double click "data" folder and select all the files in that folder

![alt text](/README_screenshots/data_select_all_ThxsPy.JPG)

Notice that alongside the data files, there's also a "sample_info.xlsx" file that looks like this

![alt text](/README_screenshots/sample_info_screenshot.JPG)

And Voila! Calculation is done and you're asked to save the output file with a file name of your choice. You can either write the ".xlsx" or not. The program will add one for you if you don't.

![alt text](/README_screenshots/save_output_ThxyPy.JPG)


## License

[BSD](https://opensource.org/licenses/BSD-2-Clause)
