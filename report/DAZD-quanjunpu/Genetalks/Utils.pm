package Genetalks::Utils;
use File::Path;

sub checkDir{
	my ($self, $dir) = @_;
	unless(-d $dir){
		mkpath($dir, 0, 0755);
	}
}





1;
