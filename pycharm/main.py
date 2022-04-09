import timeit
from collections import Counter
from collections import defaultdict

WORDS_PLASMID_TXT = r"cog_words_plasmid.txt"
WORDS_BAC_TXT = r"cog_words_bac.txt"


# take one row from DB and check if cog appears in it
def find_cog(data, cog_map, cog, l):
    for row in data:
        for i in range(5, len(row)):
            if cog == row[i]:
                find_words(cog_map, row, i, l)  # find all neighbors of cog
                break


# find all words (cog+neighbors) in a specific row
def find_words(cog_map, row, i, l):
    for l_param in range(2, l + 1):  # number of neighbours 2-10
        if l_param <= len(row) - 5:  # 5 is the index of the first cog in the row
            for j in range(5, len(row) - l_param + 1):
                word = tuple(row[j: j + l_param])  # tuple is hashable (the length of word is l_param)
                if row[i] in word:  # if the cog is in the word
                    if word in cog_map:  # if the word is already in the hashmap
                        value_of_word = cog_map[word]
                        value_of_word[0] += 1
                        organism = row[3]  # row[3] is the organism's name
                        value_of_word[1].add(organism)  # add the organism's name to the set
                        cog_map[word] = value_of_word  # set to the map the new value of key 'word'
                    else:
                        cog_map[word] = [1, {row[3]}]  # add a new key 'word' to the map


# filter the map and only consider words that appear in at least q different organisms
def bigger_than_q(cog_map, q, counter):
    for key in cog_map:
        value = cog_map[key]
        if len(value[1]) >= q:  # check if the size of the organisms set is at least q
            counter[key] = len(value[1])


def bigger_than_q2(q, counter):
    new_counter = Counter()
    for word in counter:
        count = counter[word]
        if counter[word] >= q:  # check if the size of the organisms set is at least q
            new_counter[word] = count
    return new_counter


def sort_output(counter, output):
    # Inserting the words into the output by the number of their appearances in the genome, in descending order
    for k, v in counter.most_common():
        output.append(k)

    group = defaultdict(list)
    # sort the output by the length
    for c in output:
        group[len(c)].append(c)
    return group


# part 2

def k_instances(sorted_groups, output, k, counter2):
    for word in output:  # Go over all the words in output
        len_word = len(word)
        sum_k = 1
        check_key = len_word + sum_k
        while sum_k <= k and check_key in sorted_groups.keys():  # check if the next group(by key) contain in sorted_group
            group = sorted_groups.get(check_key)
            find_word = is_find_word_k(word, group)  # number of word in each group
            counter2[word] = counter2[word] + find_word  # insert to new Counter
            sum_k += 1
            check_key += 1


def is_find_word_k(word, group):
    total_sum = 0
    for arr in group:
        if word[0] == arr[0] and word[len(word) - 1] == arr[len(arr) - 1]:  # the insertion in the beginning and in end.
            j = 1
            for i in range(1, len(arr)):
                if word[j] == arr[i]:
                    j += 1
                if j == len(word):  # so we arrive to the last word in the arr, and we check before it. so we done
                    total_sum += 1
                    break
    return total_sum


def update_counter2(counter1, counter2):
    for word in counter2:
        counter2[word] += counter1[word]

# part 3


def d_deletions(sorted_groups, output, d, counter3):
    for word in output:  # Go over all the words in output
        len_word = len(word)
        if len_word > 2:
            sum_d = 1
            check_key = len_word - sum_d
            while sum_d <= d and check_key >= 2 and check_key in sorted_groups.keys():  # check if the next group(by key) contain in sorted_group
                group = sorted_groups.get(check_key)
                find_word = is_find_word_d(word, group)  # number of word in each group
                counter3[word] = counter3[word] + find_word  # insert to new Counter
                sum_d += 1
                check_key -= 1


def is_find_word_d(word, group):
    total_sum = 0
    for arr in group:
            j = 0
            for i in range(0, len(word)):
                if word[i] == arr[j]:
                    j += 1
                if j == len(arr):  # so we arrive to the last word in the arr, and we check before it. so we done
                    total_sum += 1
                    break
    return total_sum


def update_counter3(counter2, counter3):
    for word in counter2:
        counter3[word] += counter2[word]


# open, read and split both of the files. return the words that contain the cog and appear in at least q organisms
def main(cog, q, l, k, d):
    start_time = timeit.default_timer()
    # first file - plasmid
    with open(WORDS_PLASMID_TXT, "rb") as f:
        data1 = f.readlines()  # read content as lines

    # split each line by tab (\t) and remove the new line at the last cell
    split1 = [x.split(b"\t")[:-1] for x in data1]

    # split the first cell by hashtag (#) and combine with the rest of the list
    data1 = [list([*x[0].split(b"#"), *x[1:]]) for x in split1]

    counter1 = Counter()  # create a Counter with key = word, value = the size of the group that the organism appears in it (q)
    cog_map = {}  # create a hashmap with the word (a tuple) as a key and array of [number_of_appears, set_of_organism] as the value
    find_cog(data1, cog_map, cog, l)  # a map of words that contain the cog from the first DB

    # open second file - bacteria
    with open(WORDS_BAC_TXT, "rb") as f:
        data2 = f.readlines()  # read content as lines

    # split each line by tab (\t) and remove the new line at the last cell
    split2 = [x.split(b"\t")[:-1] for x in data2]

    # split the first cell by hashtag (#) and combine with the rest of the list
    data2 = [list([*x[0].split(b"#"), *x[1:]]) for x in split2]

    find_cog(data2, cog_map, cog, l)  # a map of words that contain the cog from both DB
    bigger_than_q(cog_map, 1, counter1)  # update the counter withe words that the value bigger than q
    output = []
    sorted_groups = sort_output(counter1, output)  # a  sort list of words that appear at least one
    counter2 = Counter()
    if k == 0:  # solution from part A
        new_counter = bigger_than_q2(q, counter1)
        output = []  # to clear the output
        sorted_groups_k = sort_output(new_counter, output)
        counter2 = counter1
        # print(sorted_groups_k)
    else:
        k_instances(sorted_groups, output, k, counter2)
        update_counter2(counter1, counter2)
        new_counter = bigger_than_q2(q, counter2)
        output = []  # to clear the output
        sorted_groups_k = sort_output(new_counter, output)
        # print(sorted_groups_k)

    counter3 = Counter()
    d_deletions(sorted_groups_k, output, d, counter3)
    update_counter3(counter2, counter3)
    new_counter = bigger_than_q2(q, counter3)
    output = []  # to clear the output
    sorted_groups_d = sort_output(new_counter, output)
    print(new_counter)
    print("\n")
    print(sorted_groups_d)
    # time_stop3 = timeit.default_timer()
    # print('Time of part3: ', time_stop3-start_time)


if __name__ == '__main__':
    main(b"0586", 20, 10, 1, 1)
