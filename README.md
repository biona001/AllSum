# AllSum

This is a Julia implementation of the L0Learn algorithm proposed in [Ensembled best subset selection using summary statistics for polygenic risk](https://www.biorxiv.org/content/biorxiv/early/2023/09/27/2023.09.25.559307.full.pdf).

Given LD matrix 𝐑 (approximating the gram matrix ($𝐗^\top𝐗$) and a vector 𝐫 representing the scaled Z-score, we solve the following objective for 𝛃:

$$\frac{1}{2} 𝛃^\top𝐑𝛃 - 𝛃^\top 𝐫 + \lambda_0\|𝛃\|_0 + \lambda_1\|𝛃\|_1 + \lambda_2\|𝛃\|_2^2$$

