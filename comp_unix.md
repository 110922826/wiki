# Linux/Unix command-line tutorial
## Preface
This tutorial is based on a Linux/Unix *command-line*. Using the *command-line* requires a Linux/Unix operating system. The easiest way to try out a Linux system without actually installing it on your computer is a [LiveCD](https://en.wikipedia.org/wiki/Live_CD). A LiveCD is a DVD that you prepare (e.g. burn a Linux distribution on it) and insert in your computer. You would restart you computer and can run Linux from the DVD without any installation requirements. This is helpful for trying out a distibution of Linux not for actual work.

Another route would be to use a virtual machine. Software to create a virtual machine is free, e.g. [VirtualBox](https://www.virtualbox.org/). 

Common flavors of Linux ready for download are e.g. [Ubuntu](https://help.ubuntu.com/community/LiveCD) or if you are thinking of going the bioinformatics route, [BioLinux](http://nebc.nerc.ac.uk/tools/bio-linux/bio-linux-7-info), which includes many pre-installed bioinformatics tools.

##Introduction##
This is a collection of commands and programs I put together for working under Linux/Unix shells. It is not comprehensive. It includes very basic stuff. Tutorial style. This is bash syntax but most of it  will work on other shells (tcsh, sh) as well.

You need a editor use nedit, gedit, emacs. Editing on the shell: emacs -nw or vi.

Hint! I use the ">" sign to indicate the bash-shell prompt and "#" to indicate a
comment in code.

Open a terminal window and you are are ready to go.

**Help about a program e.g. 'ls':**

```bash
> man ls
> ls -h
```
Another very helpful resource is the [explainshell.com ](http://www.explainshell.com ) webpage,  that lets you write down a *command-line* to see the help text that matches each argument.


**Investigate directory / list directory:**

```bash
# list the current directory implicitly 
> ls
# the same in a nicer format
> ls -la
# List a particular directory (e.g. temp/) explicitly
> ls temp/
```

**Moving arround in the file system**

```bash
# what directory am I in?
> pwd

# change into directory "temp"
> cd temp/

# Go one directory up in the directory tree:
> cd ..

# Go to your home directory from any position in the directory tree:
> cd

# A shortcut for the home directory is ~/
# This command will change to /home/user/temp from any position in the directory tree
> cd ~/temp
```


## File-handling##
**Create a new empty text-file:**

```bash
> touch file1.txt
```

**Copy a file (file1.txt) to a other location or file:**

```bash
> cp file1.txt file2.txt
```

**Copy directories:**

```bash
> cp -r dir1 dir2
```

**Move a file/directory:**

```bash
> mv file1.txt file2.txt
> mv dir1 dir2
```

**Delete a file (*caution*):**

```bash
> rm file1.txt
```

**Delete a dir:**

```bash
> rm -r dir/
```

Warning! Try to avoid using "rm *". This will erase all files in the directory.

**Concatenate content of files -> will print on stdout:**

```bash
> cat file1.txt file2.txt
# all files starting with "file":
> cat file*
# print content from one file to stdout:
> cat file1.txt
```

**Look into files**

```bash
> less file1.txt
```

**Print head/tail of files**

```bash
# first 15 lines:
> head -15 file1.txt

# or:
> cat file1.txt | head -15

# last 15 lines:
> tail -15 file1.txt
> cat file1.txt | tail -15

# remove header line:
> cat file1.txt | tail -n +2 
```

Note! This is the first time we see the concept of a *pipe* "|". The pipe allows us to use the output from one program as input into another program. Many unix programs print their output by default to standard output (stdout), which is generally the shell. Most programs are also able to read their input from standard input (stdin), the shell. Thus we are able to chain commands/programs together to produce a final result. This is a very useful concept on the shell and we will use this all the time.

**Count number of rows of a file:**

```bash
> wc -l file1.txt
> cat file1.txt | wc -l
```

**Sorting files**

```bash
# sort on complete line:
> sort file1.txt
> cat file1.txt | sort

# Sort a comma-seperated file on third field:
> cat file1.txt | sort -t ',' -k 3,3

# Sort a tab-eperated file on a field (here second) but keep header line intact:
> cat file1.txt | awk 'NR==1; NR > 1 {print $0 | "sort -n -k 2,2"}' 
```

**Make lines of a file uniq / needs input from "sort":**

```bash
> cat file1.txt | sort | uniq
```

Hint! There are differnt *pipes* to direct certain outputs of a program to a
different location:

**Now we want to write the output from a chain of commands into a new file.**

```bash
# Pipe output from a program into a file:
> cat file1.txt | wc -l > new_file.txt
```

**Or append onto an existing file (if file does not exists: command also creates a file):**

```bash
> cat file1.txt | wc >> existing_file.txt
```

**Extract columns of a file:**

```bash
# less file1.txt | cut -d'seperator' -fCOLUMN,COLUMN,...   e.g.:
> cat file1.txt | cut -d ',' -f 1,2,5-7
```

**Now we will combine everything:**

1. extract from a comma seperated file the second column
2. sort the lines
3. make them uniq
4. count them
5. write number of lines into a file

```bash
> cat file1.txt | cut -d',' -f2 | sort | uniq | wc -l > number_of_lines.txt
```

**Search for pattern in a file. grep/egrep**

```bash
# print only lines of a file that contain a pattern:
> cat file1.txt | grep 'STRING'
> grep 'STRING' file1.txt

# print only lines that do _not_ contain the pattern:
> cat file1.txt | grep -v 'STRING'
```

**The same with regular expressions:**

```bash
# Here "egrep" matches every line that contains the word "REGEXP"
> cat file1.txt | egrep 'REGEXP'
# print non-matching lines:
> cat file1.txt | egrep -v 'REGEXP'
```

 
**In file substitutions using sed**

```bash
# substiute all tabs with commas (g for global)
> cat file1.txt | sed 's/\t/,/g' 
```

**Perl on command line:**

```bash
> cat file1.txt | perl -lne 'PERL-COMMAND'
```

**Transform windows newline (\r\n) into unix-newline (\n)**

```bash
> cat windows-file | tr '\015\012' '\012' > unix-file
# other way around
> cat unix-file | tr '\012' '\015\012' > windows-file
```

**Convert xls to csv**

```bash
# 1. gnumeric
> ssconvert Book1.xlsx newfile.csv

# 2. libreoffice
> libreoffice --headless --convert-to csv file.xls --outdir . 
```

##Archive and compression magic##

**First we compress a single file:**

```bash
> gzip file1.txt
# will produce a file called file1.txt.gz, and delete file1.txt
```

**We do the same on the fly reading from stdin:**

```bash
> cat file1.txt | gzip > file1.txt.gz
# Now file1.txt still exists
```

**We do not need to decompress a file to use its content (most of my text files are stored in gzip format):**

```bash
> zless file1.txt.gz
> zcat file1.txt.gz | less
> zcat file1.txt.gz | wc -l
```

**We create an archive and store two directories in it and pass it to gzip for compression**

```bash
> tar cvf - /foo /bar | gzip > foobar.tar.gz
```

**Extract  archive again:**

```bash
> tar xvzf foobar.tar.gz
```

**ADVANCED tar: Using tar to send files compressed over networks via ssh:**

```bash
# tar czf - /some/file | ssh username@computername_or_ip "tar xzf - -C /destination"
# e.g.:
> tar czf - *.pdf | ssh seb@SAUSAGE "tar xzf - -C ~/temp"
```

**USING scp: Copy a file to another computer (the colon is important)**

```bash
# scp file1.txt username@computername_or_ip:/dir_to_copy
# e.g.:
> scp file1.txt seb@SAUSAGE:~/temp

# copy whole directory to another computer
> scp -r data/ seb@SAUSAGE:~/temp

# copy from another computer over the network to your local machine
> scp seb@SAUSAGE:~/file1.txt .
```


##File and directory privileges 

**Change the privileges on directories and files.
For some of these commands to work you need root/sudo privileges.**

```bash
# Read, execute and write rights for all on all files and dirs:
> chmod 777 *

# Just for myself and the group:
> chmod 770 *

# Just for myself:
> chmod 700 *

# With 5 just read and execute rights:
> chmod 755 *

# For all files in a directory:
> chmod -R 755 dir1/*
```

##Dealing with running processes

**See processes:**

```bash
# get an overview and stats of jobs and their resource requirements 
> top
# see all system jobs
> ps -ef
# see only your jobs
> jobs
```

**Kill a process**

```bash
# kill %jobnumber/processnumber
# e.g.:
> kill %1112
# Somewhat more rude kill
> kill -9 %1112
```

**Download a file with limited download-rate**

```bash
# wget --limit-rate=12k source
> wget --limit-rate=12k http://PATH/bla.tar.gz
```

##Random stuff

**Error investigation (unix, mysql,...)**

```bash
# perror id
> perror 1234
```

**Synchronization of folders**

```bash
> rsync -avr --delete --progress --exclude "GNP" --exclude "RCS" --exclude "CVS" --exclude "CVSROOT" /home/seb/projects/ /backup_projects >> logfile.txt

# over the network to a different machine
> rsync -avr --delete --progress --exclude "GNP" --exclude "RCS" --exclude "CVS" --exclude "CVSROOT" /home/seb/projects/ seb@SAUSAGE:/backup 
```

**Start a program at specific time point**

```bash
# at -fCOMMAND HH:MM
# e.g. BASH-SCRIPT.sh holds a line like
# wget http://PATH/bla.tar.gz  
> at -fBASH-SCRIP.sh 12:45

# Look at scheduled jobs 
> atq

# Delete a scheduled job
> atrm
```

**MYSQL**

```bash
# Enter a mysq-shell:
> mysql -u seb -p

# Create new user and grant im all priviledges from any host od sanbi
mysql> create user sebastian identified by 'seb'
mysql> grant all on *.* to sebastian@"%.sanbi.ac.za" identified by 'seb';
mysql> grant all on *.* to sebastian@localhost identified by 'seb';

# Do a mysql dump on the shell
# mysqldump --databases db_name1 [db_name2 ...] > my_databases.sql
> mysqldump --databases MDB > my_databases.sql

# Load db back in mysql installation
> mysql -uroot -p < my_databases.sql
```

**Mount a USB device**

```bash
> sudo mount /dev/sda1 /mnt

# unmount
> umount /mnt
```

**_File: comp_unix.md - Sebastian Schmeier - Last update: 2014/02/14_**
