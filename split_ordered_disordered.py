####classify ordered and disordered into txt files

c= open('C:\Users\Αναστασία\Desktop\Daily_work\Raw_materials_for_samples\surf3.txt', 'r')
output1=open('C:\Users\Αναστασία\Desktop\Daily_work\Raw_materials_for_samples\ordered_text.txt', 'w+')
output2=open('C:\Users\Αναστασία\Desktop\Daily_work\Raw_materials_for_samples\disordered_text.txt', 'w+')


for line in c:
    if line[-2]=='0':
        output1.write(line)
    else:
        output2.write(line)
