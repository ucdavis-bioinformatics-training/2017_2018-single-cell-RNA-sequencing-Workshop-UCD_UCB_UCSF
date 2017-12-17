# Installing the workshop software

**ALL of this can be done on the head node, ganesh**

We are going to install the software needed for this workshop

1. FLASH2 to overlap reads (https://github.com/dstreett/FLASH2)
2. RDP to classify reads into taxon (https://github.com/rdpstaff/RDPTools)
3. dbcAmplicons pipeline for processing amplicon sequences to abundance tables (https://github.com/msettles/dbcAmplicons)

optional

4. bowtie2 (https://github.com/BenLangmead/bowtie2)

**System Requirements**

* git
* java jdk
* ant
* python virtualenv
* biom requires numpy

---

**1\.** First, create a directory for the workshop in your home directory:

    cd 
    mkdir mca_example

and two other directors 

	mkdir mca_example/src
	mkdir mca_example/bin

**2\.** Now lets add the new bin directory to our PATH in a \.bash_profile file

using your favorite text editor, _nano_ is simple, add the line

	export PATH=~/mca_example/bin:$PATH

to a file named \.bash_profile [node the leading \. as its a 'hidden' file]. Then on the command line, so that we don't have to first log out then log back in.

	source ~/.bash_profile

---

**3\.** Install **FLASH2** into src and link the exectuable into bin

	cd mca_example/src
	git clone https://github.com/dstreett/FLASH2.git
	cd FLASH2/
	make
	ln -s ~/mca_example/src/FLASH2/flash2 ~/mca_example/bin/.
	# test installation, should see help documentation
	flash2 -h
	cd ..

---	

**4\.a** Install apache ant, need for RDP

	wget http://mirrors.ibiblio.org/apache//ant/binaries/apache-ant-1.10.1-bin.tar.gz
	tar xzvf apache-ant-1.10.1-bin.tar.gz
	ln -s ~/mca_example/src/apache-ant-1.10.1/bin/ant ~/mca_example/bin/.

**4\.b** Install the **Ribsomal Database Project** (RDP) into src

	module load java/jdk1.8
	git clone https://github.com/rdpstaff/RDPTools.git
	cd RDPTools/
	git submodule init
	git submodule update
	make
	# test installation, should see help documentation for classify
	java -jar classifier.jar classify
	cd ..

**4\.c** Add the location of classifier.jar as a variable in our \.bash_profile

using your favorite text editor, _nano_ is simple, add the line

	module load java/jdk1.8
	export RDP_PATH=~/mca_example/src/RDPTools

to a file named ~/\.bash_profile, then source it	

	source ~/.bash_profile

---

**5\.a** Setup a python virtual environment for dbcAmplicons

	module load python-libs/2.7.6-ubuntu
	virtualenv dbcA_virtualenv

**5\.b** now lets set the virtual environment to activate on login by adding it to our \.bash_profile

using your favorite text editor, _nano_ is simple, add the lines

	module load python-libs/2.7.6-ubuntu
	source ~/mca_example/src/dbcA_virtualenv/bin/activate
	
to a file named ~/\.bash_profile, then source it	

	source ~/.bash_profile

---

**6\.** install **dbcAmplicons**

	pip install biom-format
	git clone https://github.com/msettles/dbcAmplicons.git
	cd dbcAmplicons/
	python setup.py install
	# test installation, should see help documentation
	dbcAmplicons -h
	cd ..

(Optional) You could also test the dbcAmplicons installation by running the script, test_dbAmplicons.sh, under the tests folds (in dbcAmplicions).
---

**Lets Review**

We created a directory for the workshop, in that directory we created two folders src and bin. We've installed FLASH2, RDP and dbcAmplicons. We've placed the executable for flash in a bin folder and added the folder to our PATH. We created an environment variable for the RDP classifier. We've created a python virtual environment and then installed the python package biom-format using pip and the dbcAmplions package using setup.py.

1. You should have verified all the software works by viewing their help docs
2. Verify that the flash executable is indeed in the bin folder
3. We've modified your PATH and added 1 new environment variable, RDP_PATH, verify the PATH and the new env variable.
4. We added a multiple lines to your \.bash_profile. How many lines?

Now log out, log back in and verify that each application still works. Ex.

	flash2 -h

To verify RDP use 

	java -jar $RDP_PATH/classifier.jar classify
	dbcAmplicons -h

You can check the current version of everything with, we usually save the output of this in the project file to remind ourselves which versions were run.

	dbcVersionReport.sh

You can test the dbcAmplicons installation buy running the script test_dbAmplicons.sh, in the tests folder.

If for some reason installation failed you can extract a copy from the workshop share

	cd
	cd mca_example
	tar xzvf /share/biocore/workshops/2017_Sept_MCA/software.tar.gz

You still need to set up the same environment variable in your \.bash_profile

---

**7\.** Last lets copy the workshop data into our home directory.

	cd 
	cd mca_example
	cp -r /share/biocore/workshops/2017_Sept_MCA/Illumina_Reads .

Take a look at the files ... what is inside the Illumina_Reads folder?

