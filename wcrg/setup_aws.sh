#!/bin/bash
#chmod u+rwx setup_aws.sh
for var in "$@"
do
	scp -i chrisstartscloud.pem -r ../eigen ubuntu@$var:/home/ubuntu/
	scp -i chrisstartscloud.pem *.cpp *.hpp makefile install.sh ubuntu@$var:/home/ubuntu/
	scp -i chrisstartscloud.pem *.pem ubuntu@$var:/home/ubuntu/.ssh/id_rsa
	ssh -i chrisstartscloud.pem ubuntu@$var -t 'sudo ./install.sh'
done

for var1 in "$@"
do
	count=0
	for var2 in "$@"
	do
  	ssh -i chrisstartscloud.pem ubuntu@$var1 -t "sudo sed -i ''1a$var2\ host$count'' /etc/hosts"
        ssh -i chrisstartscloud.pem ubuntu@$var1 -t "ssh-keyscan host$count >> /home/ubuntu/.ssh/known_hosts"
	let count++
	done
done


line_num=$(sed -n "/aws/=" /etc/ssh_config) 
line_num=`expr $line_num + 1`
sudo sed -i.bak "$line_num s/HostName.*/HostName $1/" /etc/ssh_config
