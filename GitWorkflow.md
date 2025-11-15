# Basic Git/GitHub Workflow

This repo uses a simple Git workflow. Follow these steps to work with it safely.

## 1. Clone the repository (first time only)

Open a terminal and run:

```bash
git clone https://github.com/DoggoOP/SSI_LiquidProp_PT_model.git
cd SSI_LiquidProp_PT_model
```
## 2. Always pull before you start working

Before you make any changes, sync with GitHub:
```bash
git pull
```
This gets the latest version so you don’t overwrite other people’s work.

## 3. Make your changes

Edit files, add code, etc.

## 4. Add, commit, and push your changes

When you’re ready to save your work:
```bash
git add .
git commit -m "your descriptive commit message"
git push
```

git add . stages all modified/new files.

git commit records a snapshot of your changes.

git push uploads your commits to GitHub.
