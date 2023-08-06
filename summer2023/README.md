# courses

## Notes for Jupyter

- If you're trying to start Jupyter terminal, you will be directed back to homepage. Please click the `Running` tab, and you can find the terminal you just started. Click it and you will enter the real terminal.
- If you're using SOS notebook, you may run into errors after executing the first Bash command in `Calysto Bash` block. Please wait several seconds and rerun it.

## Notes for Exercises

For Exercise of lecture IV and lecture V, we have to first copy the prepared notebook to our own Jupyter folder.

Please open jupyter, create (click `New`) a Calysto Bash Notebook, or start a terminal, and run following commands:

```bash
# Exercise for lecture IV
cp /course/popgen23/notebooks/admixExercise_popgen23.ipynb ~/

# Exercise for lecture V
cp /course/popgen23/notebooks/populationStructureII.ipynb ~/
```

After copying files, please refresh your Jupyter homepage, and you will find these two new notebooks.

## FAQ

- When I run bash command in Jupyter notebook (SOS notebook), I have this error:

    ```text
    No subkernel named Bash is found. Please use magic "%use" without option to see a list of available kernels and language modules.
    ```

    On the top right corner of current code block, you may see a frame/label called `Bash` . Please click on it and then choose the "Calysto Bash".
