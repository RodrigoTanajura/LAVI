






def peak_finder(amp, freq, intervalo=10):
    densidade = intervalo*2+1
    picos = find_peaks(amp, prominence=5e-6)
    speed = freq
    for i in range(len(picos[0])):
        speed_range = np.append(speed_range, speed[picos[0][i]-intervalo:picos[0][i]+intervalo])
        amp2 = np.append(amp2, amp[picos[0][i]-intervalo:picos[0][i]+intervalo])
    return amp2, speed_range