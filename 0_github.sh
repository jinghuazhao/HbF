
#!/usr/bin/bash

function setup()
{
  git init
  git remote add origin git@github.com:jinghuazhao/HbF.git
}

function send()
{
  git add README.md ?_* .gitignore
  git commit -m "HbF"
  git branch -M main
  git push -u origin main
}

send
