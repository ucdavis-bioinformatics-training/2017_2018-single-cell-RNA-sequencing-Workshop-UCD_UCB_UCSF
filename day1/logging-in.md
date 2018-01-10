# Logging In & Transferring Files

Nik Joshi

najoshi@ucdavis.edu


### For Macs/Linux - Logging In

1. Open a Terminal (usually under Applications/Utilities on a Mac)
2. Cut and paste this into the terminal:

        ssh [username]@ganesh.genomecenter.ucdavis.edu

   where [username] is replaced with your course username, ex. class12. Press Enter.

3. Type in your password. No characters will display when you are typing. Press Enter.

### For Macs/Linux - Transferring files

1. Using scp (secure copy, a remote file copying program):

	scp [from] [to]
	Ex.
        scp [username]@ganesh.genomecenter.ucdavis.edu:[full path to file] .

   Replace [username] with your cource username with your username and replace '[full path to file]' with the full path to the file you want to transfer. Note that there is a "." at the end of the command, which is where to put the file, i.e. your current working directory. You will have to type in your password.


### For Windows - Logging In

1. Open up PuTTY. If you haven't installed PuTTY, get it [here](http://www.putty.org/).
2. In the Host Name field, type **ganesh.genomecenter.ucdavis.edu**
3. Make sure the Connection Type is SSH.
4. Press "Open". It will ask you for your username and password.


### For Windows - Transferring files

1. Open up WinSCP. If you haven't installed it, get WinSCP [here](https://winscp.net/eng/download.php).
2. In the Host Name field, type **ganesh.genomecenter.ucdavis.edu**
3. Type in your username and password.
4. Make sure the File Protocol is SFTP.
5. Press "Login".

OR

### For MAC or Windows - Transferring files

Cyberduck - https://cyberduck.io/?l=en


