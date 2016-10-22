final = ['1', '2', '3', '4']
hmm = ['1', '2123532', '325215', '432']
# tmp = len(final)
# for a in final:
#     if a not in hmm:
#         print(a)
#     else:
#         final.remove(a)
# print(tmp, len(final))
a = [x for x in final if not determine(x)]
