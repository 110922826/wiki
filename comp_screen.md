# screen how-to

## Overview

Freely liberated from [wikipedia](https://en.wikipedia.org/wiki/): [GNU Screen](https://en.wikipedia.org/wiki/GNU_Screen) is a software application that can be used to multiplex several virtual consoles, allowing a user to access multiple separate terminal sessions inside a single terminal window or remote terminal session. It is useful for dealing with multiple programs from a command line interface, and for separating programs from the Unix shell that started the program.

**Why is this a useful tool?** In essence, whenever I use a shell, I use screen (e.g. the first command on the newly opened shell is either "screen" to start a new screen session or connect to an old screen session) and so should you! If you use screen only for one unbeatable reason, image the following scenario:

*"""Image you open a new console window and execute a script/program that needs a while to run. Now by mistake you close that console window. The program that you executed in that window on a normal shell and was running for 3 hours, dies with the closing of the window. Not so in screen. Once you start a screen session and run programs in it, you can close the window (with or without disconnecting from the screen session) and your jobs will still run. You can easily open another console and connect back to the screen session. Everything will be as you left the session (e.g. the history, etc.)."""*

##The *.screenrc* file

First of all screen can be adjusted using a "config"-file (*.screenrc**) that always gets loaded when you initiate a new screen session. 

```bash
> cd
> touch .screenrc
#open the .screenrc file in an editor
> vim .screenrc
```

Here is an example of my *.screenrc*. You can copy and paste this into your *.screenrc* if you wish so.

```
startup_message off

# lastline always on with info/date/clock
hardstatus on
hardstatus alwayslastline
hardstatus string "%{=b} Tabs: %w%=%y/%m/%d %c "

# turn visual bell on
vbell on
vbell_msg "beep"

# define a bigger scrollback, default is 100 lines
defscrollback 1024

# execute shells in screnn
screen -t s1 1
screen -t s2 2
screen -t s3 3

#screen -t tr4 4 ssh -X sebastian@tr
select 1

# Change escape character sequence to ctrl-v
defescape ^vV
escape ^vV

# detach on hangup
autodetach on

# key-bindings
#bind  prev
#bind  next

#f1 and f2, forward and back
# bindkey -k k1 prev # is help in ubuntu
bindkey -k k2 prev
bindkey -k k3 next

# makes less, vim, etc. work properly
altscreen on
```

##How-to *use* screen

```bash
# start a new screen session
> screen
# start a new named screen session
> screen -S test
# reattaching to an already running session
> screen -r [PID]
# reattaching to an already running session with name "test"
> screen -r test
# Attach to a not detached screen session.
> screen -x [PID]
# Attach to a not detached screen session with name "test"
> screen -x test
```

Make use of the control-sequence within the screen session. If you used my example *.screenrc* file the control sequence is ctrl-v otherwise the default is ctrl-a (I will use the default here). The most used commands:

```bash
# create a new shell in the screen session
ctrl-a c

# rotate through your shells
ctrl-a n

# go to a particular shell e.g. 1
ctrl-a 1

# creating a logfile of the session
ctrl-a H

# Splitting the current region in tow
ctrl-a S

# make current region the only one
ctrl-a Q

# detaching from a active screen session
ctrl-a d
```




