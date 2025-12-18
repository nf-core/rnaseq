#!/usr/bin/env bash

# Customise the terminal command prompt
echo "export PROMPT_DIRTRIM=2" >> $HOME/.bashrc
echo "export PS1='\[\e[3;36m\]\w ->\[\e[0m\\] '" >> $HOME/.bashrc
export PROMPT_DIRTRIM=2
export PS1='\[\e[3;36m\]\w ->\[\e[0m\\] '

# Update Nextflow
nextflow self-update

# Update welcome message
echo "Welcome to the nf-core/rnaseq devcontainer!" > /usr/local/etc/vscode-dev-containers/first-run-notice.txt
