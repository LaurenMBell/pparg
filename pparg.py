import model_ag
import model_ant

def main():
    print("================= PPARG DRUG EFFECTIVENESS PREDICTION =================\n")
    print("Input hypothesized drug in SMILES format.")
    target = input("(ex: COc1ccc([C@@H](O)[C@H](O)CO)cc1OC): ")
    print("\n")

    model_ant.run_model_ant(target)
    model_ag.run_model_ag(target)
    return
 

if __name__ == "__main__":
    main()