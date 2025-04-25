"""
Used this to make Figure 6 for my poster - Ben Driggs
"""

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def main():
    diet_A = [3.61, 0.64, 1.98]
    diet_C = [4.24, 0, 1.91]
    diet_F = [2.88, -0.08, 1.25]
    diet_G = [3.75, 0.06, 1.72]
    
    diets = [diet_A, diet_C, diet_F, diet_G]
    
    names = ["LP AL", "LP CR", "HP AL", "HP CR"]
    
    colors = ["blue", "red", "green", "orange"]
    shapes = ["o", "D", "*"]
    
    for i in range(len(diets)):
        for j in range(len(diet_A)):
            plt.scatter(names[i], diets[i][j], color=colors[i], marker=shapes[j], s=50)
    
    circle = Line2D([0], [0], label='Avg Essential AA Change', linestyle="", marker="o", color='black', markersize=7)
    star = Line2D([0], [0], label='Avg Combined AA Change', linestyle="", marker="*", color='black', markersize=7)
    diamond = Line2D([0], [0], label='Avg Non-Essential AA Change', linestyle="", marker="D", color='black', markersize=7)
    
    plt.ylabel("Average Change in N-Value")
    plt.legend(handles=[circle, star, diamond], bbox_to_anchor=(0.925, 1), loc='upper right')
    plt.show()
    

if __name__ == "__main__":
    main()
    