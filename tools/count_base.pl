#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;

# 解析命令行参数
my $seperate = 0;
GetOptions('s' => \$seperate) or die "Usage: $0 <input.fasta> [-s]\n";

# 检查输入文件
my $input_file = $ARGV[0] or die "Usage: $0 <input.fasta> [-s]\n";
open my $file, "<", $input_file or die "Cannot open file: $!\n";

# 初始化变量
my ($length, $N_length, $seq_count) = (0, 0, 0);
my ($id, $seq);

# 逐行读取文件
while (<$file>) {
    chomp;
    if (/^>/) {  # 如果是 header 行
        # 处理前一个序列
        if (defined $id && defined $seq) {
            process_sequence($id, $seq, \$length, \$N_length, \$seq_count, $seperate);
        }
        # 开始新的序列
        ($id = $_) =~ s/^>//;  # 去掉开头的 ">"
        $id = (split)[0];      # 去掉可能的注释部分
        $seq = "";             # 重置序列
    } else {
        # 拼接序列
        $seq .= $_;
    }
}

# 处理最后一个序列
if (defined $id && defined $seq) {
    process_sequence($id, $seq, \$length, \$N_length, \$seq_count, $seperate);
}

# 输出总结果
my $missing = $N_length / $length;
print "All\t$length\t$N_length\t$missing\t$seq_count\n";

# 关闭文件
close $file;

# 处理单个序列的子程序
sub process_sequence {
    my ($id, $seq, $length_ref, $N_length_ref, $seq_count_ref, $seperate) = @_;
    $seq =~ s/\s+//g;  # 删除空白字符
    my $N_count = $seq =~ tr/Nn//;
    my $chr_len = length $seq;
    $$length_ref += $chr_len;
    return unless $chr_len > 0;
    $$seq_count_ref++;
    my $chr_mis = $N_count / $chr_len;
    print "$id\t$chr_len\t$N_count\t$chr_mis\n" if $seperate;
    $$N_length_ref += $N_count;
}
