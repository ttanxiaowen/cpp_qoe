import random
facility = ["shanghai", "nanjing", "beijing", "hainan", "lasa", "ulu" ]
choices_f = [1/6]*6
kbps = [2,4,8,10,15,18,20,25,30,50]
ratio= [0.1]*10
choices_K = random.choices(kbps,weights=ratio,k=2000)
choices_F = random.choices(facility,weights=choices_f,k=2000)
print("good")
with open("traffic_new.txt","w") as file:
    for i in range(2000):
        random_choice = random.choice([0,1])
        if random_choice == 1:
            file.write(f"'src_name':car{i},'target_name':{choices_F[i]},'needrate':{choices_K[i]*1e6}\n")
        else:
            file.write(f"'src_name':{choices_F[i]},'target_name':car{i},'needrate':{choices_K[i]*1e6}\n")