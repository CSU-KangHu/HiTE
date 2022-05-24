# Definition for a binary tree node.
class TreeNode:
    def __init__(self, x):
        self.val = x
        self.left = None
        self.right = None


def numIslands(grid):
    row = len(grid)  # 行数
    col = len(grid[0])  # 列数
    for i in range(row):
        for j in range(col):
            path = []
            dfs(grid, i, j, row, col, path)
            print(path)

t=15

def dfs(grid, x, y, row, col, path):
    cur_node = grid[x][y]
    #print(cur_node)
    path.append((x, cur_node))
    # 递归的终止条件
    if x < 0 or x >= row-1:
        return
    # all cell in next row
    for j in range(len(grid[x+1])):
        child = grid[x+1][j]
        if child-cur_node > 0 and child-cur_node < t:
            dfs(grid, x+1, j, row, col, path)

if __name__ == '__main__':
    grid = [[0, 10, 20, 30], [10, 40, 80, 100], [5, 15, 50, 70], [20, 35, 55, 60]]
    numIslands(grid)
