def divide_sequences(text_sequence, nr, ninr):
    out_text_list = []
    counter = ninr
    for i in range(0,len(text_sequence), nr):
        subtexto = text_sequence[i:i+nr]
        tmp_list =[]
        for j in range(0,len(subtexto), ninr):
            tmp_list.append(subtexto[j:j+ninr])
        text_out = ' '.join(tmp_list)
        l_index = []
        for i in range(int(nr/ninr)):
            if ((counter + (ninr * i)) <= (len(text_sequence)+(ninr-1))):
                l_index.append(str(counter + (ninr * i)).rjust(10))
        out_text_list.append(' '.join(l_index))
        out_text_list.append(text_out)
        counter += nr
    return(out_text_list)
