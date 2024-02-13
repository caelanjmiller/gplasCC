from rich import print as rprint

def verbose_print(message, silent, status):
    if silent:
        return
    if status == 'normal':
        return print(message)
    elif status == 'start':
        return rprint(f"[bold yellow]{message}...[/bold yellow]", end="\r")
    elif status == 'end':
        return rprint(f"[bold green]{message} completed![/bold green]")
    elif status == 'error':
        rprint(f"[bold red]{message}[/bold red]")