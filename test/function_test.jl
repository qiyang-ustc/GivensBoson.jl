import GivensBoson.find_offdiagonal_maximun
@testset "Functional Test" begin
    value,index,flag = find_offdiagonal_maximun([1 0.9;0.9 1])
    @test (value,index,flag)==(0.9,[1,2],1)

    value,index,flag = find_offdiagonal_maximun([1 1.0;1.0 1])
    @test (value,index,flag)==(1.0,[1,2],1)


    value,index,flag = find_offdiagonal_maximun([1 1.0;1.0 1],zeromode=true)
    @test (value,index,flag)==(1.0,[1,2],2)

    value,index,flag = find_offdiagonal_maximun([1 .0;0 1])
    @test (value,index,flag)==(0.0,[1,1],3)
end