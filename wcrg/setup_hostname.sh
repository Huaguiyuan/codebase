#!/bin/bash
#chmod u+rwx setup_hostname.sh
line_num=$(sed -n "/aws/=" /Users/chplatt/Documents/test_ssh.txt) 
line_num=`expr $line_num + 1`
sudo sed -i.bak "$line_num s/HostName.*/HostName $1/" /etc/ssh_config
