function printRK4Steps( k1, Yk2, k2, Yk3, k3, Yk4, k4, Yfinal )

    fprintf("\nk1: [ %d, %d, %d ] norm: %d \n",...
        k1(11), k1(12), k1(13), norm(k1(11:13)))
    fprintf("Yk2: [ %d, %d, %d ] norm: %d \n",...
        Yk2(11), Yk2(12), Yk2(13), norm(Yk2(11:13)))
    fprintf("k2: [ %d, %d, %d ] norm: %d \n",...
        k2(11), k2(12), k2(13), norm(k2(11:13)))
    fprintf("Yk3: [ %d, %d, %d ] norm: %d \n",...
        Yk3(11), Yk3(12), Yk3(13), norm(Yk3(11:13)))
    fprintf("k3: [ %d, %d, %d ] norm: %d \n",...
        k3(11), k3(12), k3(13), norm(k3(11:13)))
    fprintf("Yk4: [ %d, %d, %d ] norm: %d \n",...
        Yk4(11), Yk4(12), Yk4(13), norm(Yk4(11:13)))
    fprintf("k4: [ %d, %d, %d ] norm: %d \n",...
        k4(11), k4(12), k4(13), norm(k4(11:13)))
    fprintf("Yfinal: [ %d, %d, %d ] norm: %d \n",...
        Yfinal(11), Yfinal(12), Yfinal(13), norm(Yfinal(11:13)))

end
