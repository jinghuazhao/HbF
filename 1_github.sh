
#!/usr/bin/bash

echo "# HbF" >> README.md
git init
git add README.md ?_*
git commit -m "first commit"
git branch -M main
git remote add origin git@github.com:jinghuazhao/HbF.git
git push -u origin main
