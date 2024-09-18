# OPACOS

The framework developed is useful in order to get the probability distribution of parameter of stars in Milky Way as observed by us synthetically for a section of our Sky. Using above framework we can get the ditribution function of parameters for a section of sky which can be used as a prior in Galactic surveys.

Our framework can be used to get the probability distribution of metallicity or age given priors of (age,distance) or (metallicity, distance) respectively. It accounts for selection effects and generates a probabilty distribution of parameters as would be seen by us for a section of sky.

The above framework is divided into two components:
- The first component which uses the extinction code ALextin and the probability calculation code that we created and calculates the probability of observing stars.
- The second component which uses the probability along with the GALAXIA code results in order to simulate the parameter distribution of stars that we would observe. Here GALAXIA is a simulation code which simulates the Galaxy for us which we sample for pur use.

<h2>Installation:</h2>

Download the zip file in the desired folder.

- Unzip the file using command:

   ```tar gz -xvf OPACOS-main.zip```
- After doing it follow the documentation as mentioned ahead in-order to install the GALAXIA code: https://galaxia.sourceforge.net/
- Open the **galaxy1.py** file and change the data variable to:

   ```data = ebf.read('/user/GalaxiaData/Examples/galaxy1.ebf','/')```
   
  here instead of user you will set the location of the GalaxiaData folder that you installed in the last step. The **galaxy1.py** file uses the GALAXIA code results stored in **galaxy1.ebf** files , reads it and stores the appropriate parameters useful for us in the **galaxia.dat** file.
  
- Refer to **Survey_DATA.dat** file and install the data files as mentioned in it as it would be used to extract data to be worked on. Install both the files in the main directory that is just outside **Survey_DATA.dat** file.

<h2>Running OPACOS:</h2>

To Run OPACOS it is quite straightforward, although sooner we will further streamline the process. To run OPACOS we need to do folow following steps:

- Generate the stars in Milky Way galaxy using GALAXIA. To do that go inside GalaxiaData folder that must have been created as mention in GALAXIA doxumentation. Then go to Examples and open myparameterfile in it. Once it is open set the location of the center of the cone that you want to simulate the stars in along with the area that you are interested in. 
- Make sure that the circular patch option is selected instead of all sky survey in myparameter file.
- After setting that run ```galaxia -r myparameterfile``` command in terminal inside the Examples folder to generate stars of galaxy which is stored in galaxy1.ebf file.
- Now you need to  

 
# OPACOS
