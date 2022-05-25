# Definition for a binary tree node.
import codecs
import json
from Util import read_fasta, store_fasta

class TreeNode:
    def __init__(self, x):
        self.val = x
        self.left = None
        self.right = None


def getLongestPath(grid, row, visited, skip_threshold, region_path):
    #longest_path = {}
    for i in range(row):
        for j in range(len(grid[i])):
            # start from each point in Matrix, try dfs to find a valid path
            path = {}
            dfs(grid, i, j, row, path, visited, skip_threshold)
            #print(path)
            region_path.append(path)
    #         if len(path) > len(longest_path):
    #             longest_path = path
    # print('longest path = ' + str(longest_path))


def dfs(grid, x, y, row, path, visited, skip_threshold):
    # if current node visited, return
    if visited[x][y]:
        return
    cur_node = grid[x][y]
    # add current node to path
    # if not path.__contains__(x):
    #     path[x] = []
    # p = path[x]
    # p.append(cur_node)
    path[x] = cur_node
    visited[x][y] = True
    # current node reach to the last row
    if x >= row-1:
        return
    # all cell in next row
    for j in range(len(grid[x+1])):
        # Pruning, nodes that have been visited do not need to be accessed repeatedly
        if visited[x+1][j]:
            continue
        child_node = grid[x+1][j]
        # child node is closed to current node, then there is a path between current node and child node
        if (child_node[0] - cur_node[1]) > 0 and (child_node[0] - cur_node[1]) < skip_threshold:
            dfs(grid, x+1, j, row, path, visited, skip_threshold)

if __name__ == '__main__':
    # get matrix
    reference = '/public/home/hpc194701009/Ref/dmel-all-chromosome-r5.43.fasta'
    tmp_output_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/dmel/CRD.2022-05-24.9-36-59'
    pathMatrix_file = tmp_output_dir + '/pathMatrix.csv'
    file = open(pathMatrix_file, 'r')
    js = file.read()
    pathMatrix = json.loads(js)
    skip_threshold = 200

    # Step1: According to fragment position matrix, get all valid paths
    region_paths = {}
    # go through each Matrix, compute the valid path
    for region_index in pathMatrix:
        cur_region_matrixs = pathMatrix[region_index]
        if not region_paths.__contains__(region_index):
            region_paths[region_index] = []
        region_path = region_paths[region_index]
        for ref in cur_region_matrixs.keys():
            cur_ref_matrix = cur_region_matrixs[ref]
            row = len(cur_ref_matrix)
            # define visited matrix, avoid duplicate visits
            visited = [[False for j in range(len(cur_ref_matrix[i]))] for i in range(row)]
            #print('region id=%s, ref=%s' %(region_index, ref))
            # core
            getLongestPath(cur_ref_matrix, row, visited, skip_threshold, region_path)
        region_path.sort(key=lambda x: (len(x)))
        region_paths[region_index] = region_path

    # store region_paths for testing
    region_paths_file = tmp_output_dir + '/region_paths.csv'
    with codecs.open(region_paths_file, 'w', encoding='utf-8') as f:
        json.dump(region_paths, f)

    connected_regions_file = tmp_output_dir + '/connected_regions.csv'
    file = open(connected_regions_file, 'r')
    js = file.read()
    connected_regions = json.loads(js)

    # Step2: record all fragments in each region, including start and end position
    # generate region fragments
    # region_fragments = {region_id: {0: (f1,s1,e1), 1: (f2,s2,e2)}}
    region_fragments = {}
    # region_fragments1 = {region_id: {f1: (s1,e1), f2: (s2,e2)}}
    region_fragments1 = {}
    for ref in connected_regions.keys():
        region_dict = connected_regions[ref]
        for region_index in region_dict.keys():
            if not region_fragments.__contains__(region_index):
                region_fragments[region_index] = {}
            region_fragment = region_fragments[region_index]

            if not region_fragments1.__contains__(region_index):
                region_fragments1[region_index] = {}
            region_fragment1 = region_fragments1[region_index]

            region_frags = region_dict[region_index]
            for index, frag_item in enumerate(region_frags):
                region_fragment[index] = frag_item

            for index, frag_item in enumerate(region_frags):
                frag_name = frag_item[0]
                region_fragment1[frag_name] = (frag_item[1], frag_item[2])

    # store region_fragments for testing
    region_fragments_file = tmp_output_dir + '/region_fragments.csv'
    with codecs.open(region_fragments_file, 'w', encoding='utf-8') as f:
        json.dump(region_fragments, f)

    # store region_fragments1 for testing
    region_fragments1_file = tmp_output_dir + '/region_fragments1.csv'
    with codecs.open(region_fragments1_file, 'w', encoding='utf-8') as f:
        json.dump(region_fragments1, f)


    # Step3: According to valid paths, the fragments in the region are connected across gaps.
    # Fragments in region are connected according to the longer path in valid paths.

    # keeped_path = {region_id: [f1f2f3, f4, f5f6]}
    keeped_paths = {}
    for region_index in region_paths.keys():
        region_fragment = region_fragments[region_index]
        for path in region_paths[region_index]:
            isPathExist = True
            pathName = ''
            for i, frag_index in enumerate(path.keys()):
                if region_fragment.__contains__(frag_index):
                    frag_item = region_fragment[frag_index]
                    frag_name = frag_item[0]
                    if i == 0:
                        pathName += frag_name
                    else:
                        pathName += ',' + frag_name
                else:
                    isPathExist = False
                    break

            if isPathExist:
                if not keeped_paths.__contains__(region_index):
                    keeped_paths[region_index] = []
                paths = keeped_paths[region_index]
                paths.append(pathName)
                for frag_index in path.keys():
                    del region_fragment[frag_index]


    # store keeped_paths for testing
    keeped_paths_file = tmp_output_dir + '/keeped_paths.csv'
    with codecs.open(keeped_paths_file, 'w', encoding='utf-8') as f:
        json.dump(keeped_paths, f)


    # connected_frags = {region_id: [(f1,s1,e1), (f2f3f4,s2,e2), (f5,s3,e3)]}
    connected_frags = {}
    for region_index in keeped_paths.keys():
        keeped_frag_name = []
        paths = keeped_paths[region_index]
        region_fragment1 = region_fragments1[region_index]
        for path in paths:
            if len(path) <= 0:
                continue
            # connect path
            last_connected_frag_start = -1
            last_connected_frag_end = -1
            for i, frag_name in enumerate(path.split(',')):
                keeped_frag_name.append(frag_name)
                frag_item = region_fragment1[frag_name]
                if last_connected_frag_start == -1:
                    last_connected_frag_start = frag_item[0]
                last_connected_frag_end = frag_item[1]

            if not connected_frags.__contains__(region_index):
                connected_frags[region_index] = []
            connected_frag = connected_frags[region_index]
            connected_frag.append((path, last_connected_frag_start, last_connected_frag_end))

        for frag_name in region_fragments1[region_index]:
            frag_item = region_fragments1[region_index][frag_name]
            if frag_name not in keeped_frag_name:
                if not connected_frags.__contains__(region_index):
                    connected_frags[region_index] = []
                connected_frag = connected_frags[region_index]
                connected_frag.append((frag_name, frag_item[0], frag_item[1]))

        connected_frag = connected_frags[region_index]
        connected_frag.sort(key=lambda x: (x[1], x[2]))

    # store connected_frags for testing
    connected_frags_file = tmp_output_dir + '/connected_frags.csv'
    with codecs.open(connected_frags_file, 'w', encoding='utf-8') as f:
        json.dump(connected_frags, f)

    connected_frags_file = tmp_output_dir + '/connected_frags.csv'
    file = open(connected_frags_file, 'r')
    js = file.read()
    connected_frags = json.loads(js)

    print('total region: %d' %len(connected_frags))
    connect_region_num = 0
    for region_index in connected_frags.keys():
        for connected_frag in connected_frags[region_index]:
            name_parts = connected_frag[0].split(',')
            if len(name_parts) > 1:
                connect_region_num += 1
                break
    print('connected region: %d' % connect_region_num)

    refNames, refContigs = read_fasta(reference)
    repeats_connected_file = tmp_output_dir + '/repeats_connected.fa'
    repeats_connected = {}
    index = 0
    for region_index in connected_frags.keys():
        for connected_frag in connected_frags[region_index]:
            query_name = 'R' + str(index)
            frag_name = connected_frag[0].split(',')[0]
            ref_name = frag_name.split('-s_')[1].split('-')[0]
            seq = refContigs[ref_name][connected_frag[1]: connected_frag[2]+1]
            index += 1
            repeats_connected[query_name] = seq
    store_fasta(repeats_connected, repeats_connected_file)


